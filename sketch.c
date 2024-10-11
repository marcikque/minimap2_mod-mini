#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"
#include "math.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

typedef struct {		// kmer construction information, constant
	int k;
	const char *str;	// sequence
	int len;			// sequence length
	uint32_t rid;		// reference ID
	int is_hpc;			// enable homopolymer-compression
} kmer_info;

typedef struct {		// hpc handling information
	int is_hpc;			// propagates hpc enabling
	mm128_t *buf;		// pointer to ring buffer
	int *buf_pos;		// pointer to variable holding current pos in ring buffer
	int t;				
	int *min_pos;		// pointer to variable holding pos of current minimizer within the buffer
} hpc_info;

int static mm_calc_idx(int w, int i, int k, mm128_t min, int w_start, hpc_info *hpc_vars);

void mm_push_kmer(void *km, mm128_v *p, int glbl_index, kmer_info *kmer_vars);

/**
 * Find (w,k)-minimizers on a DNA sequence using the mod-minimizer scheme:
 * 		"The mod-minimizer: a simple and efficient sampling algorithm for long k-mers"
 * 		Ragnar Groot Koerkamp, Giulio Ermanno Pibiri
 *		bioRxiv 2024.05.25.595898; doi: https://doi.org/10.1101/2024.05.25.595898
 * 
 * Notation:
 *  tmer : newly constructed tmer by removing the first base and appending the new base
 *  W : set of tmers in current window
 *  tmer_i : i-th tmer of the window
 *  kmer_i : k-th kmer of the window
 *  h(a) : hash of a
 *  rc(a) : reverse complement of a
 *  pos_W(a) : a's position within window W (0-indexed)
 *  M : list of minimizers
 *
 * Procedure: using a sliding window, construct the entering tmer
 * 	1. W = W \ {W_0}
 * 	2. def. info = min( h(tmer) , h(rc(tmer)) )
 * 	3. W = W + {info}
 * 	4. def. min = min(W)
 * 	5. def. p = pos_W(min) mod w
 * 	6. def. kmer = min( h(kmer_p)) , h(rc(kmer_p)) )
 * 	7. M = M + {kmer}
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	// k >= r should hold
	int r = 4; //(int) ceil(log2(w+k-1) / 2) + 1; //dynamic way with needed threshold
	int t = r + ((k-r) % w);

	uint64_t shift1 = 2 * (t - 1);
	uint64_t mask = (1ULL<<2*t) - 1;
	uint64_t tmer[2] = {0,0};
	int l = 0; // bases processed since last reset/start
	int w_start = 0; // starting index of the current window
	int buf_pos = 0;
	int min_pos = 0;
	int tmer_span = 0; // bases included in the current tmer (default t, changes under HPC)
	int prev_push = -1;
	mm128_t buf[256];
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	// creating the necessary structs to call mm_push_kmer to improve readability and reduce duplicate code
	hpc_info hpc_vars = {
		.is_hpc = is_hpc,
		.buf = buf,
		.buf_pos = &buf_pos,	
		.t = t,
		.min_pos = &min_pos
	};
	kmer_info kmer_vars = {
		.k = k,
		.str = str,
		.len = len,
		.rid = rid,
		.is_hpc = is_hpc
	};

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56bits for kmers, 8bits for metadata

	memset(buf, 0xff, (w + k - t) * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (int i = 0; i < len; ++i) {
		mm128_t info = { UINT64_MAX, UINT64_MAX }; 
		int c = seq_nt4_table[(uint8_t)str[i]];

		// keeps track of the start of the current window, not skewed by HPC
		if(l >= w + k - 1) {
			++w_start;
		}

		if (c < 4) { // not an ambiguous base
			if (is_hpc) {
				// next base same? -> loop till different
				// e.g. GTAAAACTG -> GTACTG
				int skip_len = 1;
				// case: next base is identical to current
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len) {
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c) {
							break;
						}
					}
					i += skip_len - 1;
				}
				tq_push(&tq, skip_len);
				tmer_span += skip_len; // count number of bases within tmer, can be >t
				if (tq.count > t) {
					tmer_span -= tq_shift(&tq);
				} 
			} else {
				tmer_span = l + 1 < t? l + 1 : t;
			}

			tmer[0] = (tmer[0] << 2 | c) & mask;           // forward t-mer: remove first base, append new base
			tmer[1] = (tmer[1] >> 2) | (3ULL^c) << shift1; // reverse t-mer: remove last base, prepend complement of new base (A<->T, C<->G)

			if (tmer[0] == tmer[1]) {
				continue; // symmetric tmer, strand unknown
			}
			int z = tmer[0] < tmer[1]? 0 : 1; // strand: forward or backward (lexicographically smaller is chosen)
			++l;

			if (l >= t && tmer_span < 256) {
				info.x = hash64(tmer[z], mask) << 8 | tmer_span; 
				// tmer.x[63:8] = tmer[z] hashed, trimmed
				// tmer.x[7:0]  = tmer_span
				info.y = i - tmer_span + 1;
				// tmer.y[63:0] = global starting position
			}
		
		// resetting, essentially splitting the sequence, at unambigious bases
		} else {
			l = 0;
			tq.count = 0;
			tq.front = 0;
			tmer_span = 0;
			buf_pos = 0;
			min_pos = 0;
			min = (mm128_t){ UINT64_MAX, UINT64_MAX };
		}

		buf[buf_pos] = info;

		// case: new minimum found
		if (info.x < min.x) {
			min = info;
			min_pos = buf_pos;
			
			if (l >= w + k) { // check if we fully processed the first window
				prev_push = mm_calc_idx(w, i, k, min, w_start, &hpc_vars);
				mm_push_kmer(km, p, prev_push, &kmer_vars);
			}

		// case: old minimum left the window 
		} else if (buf_pos == min_pos && l >= w + k) { // && check if we fully processed the first window
			min.x = UINT64_MAX;

			// scan for the leftmost minimal element
			int w_relative;
			for (w_relative = 0; w_relative < w + k - t; ++w_relative) {
				int wrapped_index = (buf_pos + 1 + w_relative) % (w + k - t);

				if (buf[wrapped_index].x < min.x) {
					min = buf[wrapped_index];
					min_pos = wrapped_index;
				}
			}

			prev_push = mm_calc_idx(w, i, k, min, w_start, &hpc_vars);
			mm_push_kmer(km, p, prev_push, &kmer_vars);
			
		// case: consecutive windows share minimizing tmer
		} else {
			int next_push = mm_calc_idx(w, i, k, min, w_start, &hpc_vars);

			// check if the tmer mod-samples a new kmer 
			if (prev_push != next_push && l >= w + k) { // && check if we fully processed the first window
				prev_push = next_push;
				mm_push_kmer(km, p, prev_push, &kmer_vars);
			}
		}

		// case: 1st window is fully processed, not partial anymore, we need to push its minimizer
		if (l == w + k - 1) {
			prev_push = mm_calc_idx(w, i, k, min, w_start, &hpc_vars);
			mm_push_kmer(km, p, prev_push, &kmer_vars);
		}

		//inc buffer
		if (++buf_pos == w + k - t) {
			buf_pos = 0;
		}
	}
}


/** 
 * Calculates the starting index of the kmer to push to the list of minimizers
 * Scheme: mod-minimizer, minimizing tmer at pos x, window starts at pos s -> push kmer at s + ((x - s) mod w)
 * Note: HPC makes the normal scheme unusable as compression skews sequence indices
 * @param w 		window length
 * @param i 		current sequence processing index
 * @param k 		kmer length
 * @param min 		minimizing tmer in the current window
 * @param w_start 	index of the start of the window
 * @param hpc_vars 	contains information to support HPC
*/
int static mm_calc_idx(int w, int i, int k, mm128_t min, int w_start, hpc_info *hpc_vars) {
	// case: hpc is disabled
	if (!hpc_vars->is_hpc) {
		// the kmer to sample is at position: w_start + ((tmer_pos - w_start) mod w)
		int glbl_index = w_start + ((min.y - w_start) % w);
		return glbl_index;
	}

	// case: hpc enabled
	// problem: 	window starting index cannot be calculated via formula as consecutive identical bases are ignored
	// solution: 	window starting index calculation by incrementing w_start instead of using a formula
	// problem: 	position ws + ((p - ws) mod w) doesnt correspond to the right kmer anymore as compression skews the indices (j-th kmer doesnt have to start at index j-1)
	// solution: 	tmer's within ring buffer contain correct starting position and numbering (j-th tmer starts at index j-1) -> do buffer calculations

	// 1. calculate the minimizer's position within the current window
	// we know that buf_pos points to the last tmer within the current window
	int window_pos;
	if (*hpc_vars->min_pos > *hpc_vars->buf_pos) {
		window_pos = *hpc_vars->min_pos - *hpc_vars->buf_pos - 1; 
	} else {
		window_pos = (w + k - hpc_vars->t) - *hpc_vars->buf_pos + *hpc_vars->min_pos - 1;
	}

	// 2. modulo according to our scheme
	window_pos %= w;

	// 3. calculate the index within the buffer of the tmer at window_pos in our current window
	// buf[buf_pos + 1] always contains the very first tmer of our window
	int target_pos = (*hpc_vars->buf_pos + 1 + window_pos) % (w + k - hpc_vars->t);

	// 4. retrieve corresponding tmer
	mm128_t target = hpc_vars->buf[target_pos];

	// 5. retrieve the starting index
	return target.y;
}


/**
 * Constructs and pushes the kmer starting at the given index
 *
 * @param km 			thread-local memory pool; using NULL falls back to malloc()
 * @param p 			list of minimizers
 * @param glbl_index 	index of kmer to push
 * @param kmer_vars 	information needed for kmer construction
 */
void mm_push_kmer(void *km, mm128_v *p, int glbl_index, kmer_info *kmer_vars)
{
	assert(glbl_index >= 0 && glbl_index < kmer_vars->len); // should never happen

	uint64_t shift1 = 2 * (kmer_vars->k - 1);
	uint64_t mask = (1ULL<<2*kmer_vars->k) - 1; 
	uint64_t kmer[2] = {0,0};
	int kmer_span = 0; // bases included in the current kmer (default k, changes under HPC)
	int l = 0;
	int j;

	for (j=glbl_index; j<kmer_vars->len && l<kmer_vars->k; ++j) {
		int c = seq_nt4_table[(uint8_t)kmer_vars->str[j]];

		if (c>=4) { // ambiguous base
			return;
		}

		if (kmer_vars->is_hpc) {
			int skip_len = 1;
			if (j + 1 < kmer_vars->len && seq_nt4_table[(uint8_t)kmer_vars->str[j + 1]] == c) { 

				for (skip_len = 2; j + skip_len < kmer_vars->len; ++skip_len) {
					if (seq_nt4_table[(uint8_t)kmer_vars->str[j + skip_len]] != c) {
						break;
					}
				}
				
			}
			kmer_span += skip_len;
		} else {
			++kmer_span;
		}

		kmer[0] = (kmer[0] << 2 | c) & mask;       		// forward kmer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1;  // reverse tmer: reverse & A<->T, C<->G
		++l;
	}

	if (j == kmer_vars->len) { // loop ended because of j<len condition: edge case if kmer ends at last pos in sequence
		--j;
	}

	if (kmer[0] == kmer[1]) { // symmetric kmer, strand unknown
		return;
	}

	if (l >= kmer_vars->k && kmer_span < 256) {
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z = kmer[0] < kmer[1]? 0 : 1;

		info.x = hash64(kmer[z], mask) << 8 | kmer_span; 
		// kmer.x[63:8] = kmer[z] hashed, trimmed
		// kmer.x[7:0]  = kmer_span
		info.y = (uint64_t)kmer_vars->rid<<32 | (uint32_t)j<<1 | z;
		// kmer.y[63:32] = sequence ref. ID
		// kmer.y[31:1]  = last position
		// kmer.y[0]     = strand
		assert(j >= 0 && j < kmer_vars->len); // replaces asserts in assign.c (more comprehensive)
		kv_push(mm128_t, km, *p, info);
	}
}