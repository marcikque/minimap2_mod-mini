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

typedef struct {		// kmer construction information
	int k;
	const char *str;	// sequence
	int len;			// sequence length
	uint32_t rid;		// reference ID
	int is_hpc;			// enable homopolymer-compression
} kmer_info;

int static mm_calc_idx(int w, int i, int k, mm128_t min);

void mm_push_kmer(void *km, mm128_v *p, int glbl_index, kmer_info *kmer_vars);

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
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
	// "The mod-minimizer: a simple and efficientsampling algorithm for long k-mers"
	int r = 4; //(int) ceil(log2(w+k-1) / 2) + 1; // r=4 default | TODO: Threshold experiment
	int t = r + ((k-r) % w); // Section 4.2 (11:13), Def. 16
	assert(k >= r); 

	uint64_t shift1 = 2 * (t - 1);
	uint64_t mask = (1ULL<<2*t) - 1;
	uint64_t tmer[2] = {0,0};
	int l = 0;
	int buf_pos = 0;
	int min_pos = 0;
	int tmer_span = 0;
	int prev_push = -1;
	mm128_t buf[256];
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56bits for kmers, 8bits for metadata

	memset(buf, 0xff, (w + k - t) * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (int i = 0; i < len; ++i) {
		mm128_t info = { UINT64_MAX, UINT64_MAX }; 
		int c = seq_nt4_table[(uint8_t)str[i]];

		if (c < 4) { // not an ambiguous base
			if (is_hpc) {
				// next base samne? -> loop till different
				// e.g. GTAAAACTG -> GTACTG
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len) {
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c) {
							break;
						}
					}
					i += skip_len - 1;
				}
				tq_push(&tq, skip_len);
				tmer_span += skip_len;
				if (tq.count > t) {
					tmer_span -= tq_shift(&tq);
				} 
			} else {
				tmer_span = l + 1 < t? l + 1 : t;
			}

			tmer[0] = (tmer[0] << 2 | c) & mask;           // forward t-mer
			tmer[1] = (tmer[1] >> 2) | (3ULL^c) << shift1; // reverse t-mer: reverse & A<->T, C<->G

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
		} else {
			l = 0;
			tq.count = 0;
			tq.front = 0;
			tmer_span = 0;
			min = (mm128_t){ UINT64_MAX, UINT64_MAX };
			buf_pos = 0;
			min_pos = 0;
		}

		buf[buf_pos] = info;

		// new minimum found
		if (info.x < min.x) {
			min = info;
			min_pos = buf_pos;
			// check if we fully processed the first window
			if (l >= w + k) {
				prev_push = mm_calc_idx(w, i, k, min);
				kmer_info kmer_vars = {
						.k = k,
						.str = str,
						.len = len,
						.rid = rid,
						.is_hpc = is_hpc
					};
				mm_push_kmer(km, p, prev_push, &kmer_vars);
			}

		// old minimum left the window &&
		// check if we fully processed the first window
		} else if (buf_pos == min_pos && l >= w + k) {
			min.x = UINT64_MAX;

			// scan for the leftmost minimal element
			int window_relative;
			for (window_relative = 0; window_relative < w + k - t; ++window_relative) {
				int wrapped_index = (buf_pos + 1 + window_relative) % (w + k - t);

				if (buf[wrapped_index].x < min.x) {
					min = buf[wrapped_index];
					min_pos = wrapped_index;
				}
			}

			prev_push = mm_calc_idx(w, i, k, min);
			kmer_info kmer_vars = {
					.k = k,
					.str = str,
					.len = len,
					.rid = rid,
					.is_hpc = is_hpc
				};
			mm_push_kmer(km, p, prev_push, &kmer_vars);
			
		// consecutive windows have the same minimizing tmer
		} else {
			int next_push = mm_calc_idx(w, i, k, min);

			// check if the tmer mod-samples a new kmer &&
			// check if we fully processed the first window
			if (prev_push != next_push && l >= w + k) { 
				prev_push = next_push;
				kmer_info kmer_vars = {
						.k = k,
						.str = str,
						.len = len,
						.rid = rid,
						.is_hpc = is_hpc
					};
				mm_push_kmer(km, p, prev_push, &kmer_vars);
			}
		}

		//1st Window special case: push leftmost minimizer as soon as whole window is processd
		if (l == w + k - 1) {
			prev_push = mm_calc_idx(w, i, k, min);
			kmer_info kmer_vars = {
					.k = k,
					.str = str,
					.len = len,
					.rid = rid,
					.is_hpc = is_hpc
				};
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
 * @param w window length
 * @param i current sequence processing index
 * @param k kmer length
 * @param min minimizing tmer in the current window
*/
// FIXME: make HPC ready
int static mm_calc_idx(int w, int i, int k, mm128_t min) {
	int w_start = i - (w + k - 2);
	int glbl_index = w_start + ((min.y - w_start) % w);

	return glbl_index;
}


/**
 * Constructs and pushes the kmer starting at the given index
 *
 * @param km 			thread-local memory pool; using NULL falls back to malloc()
 * @param p 			list of minimizers
 * @param glbl_index 	index of kmer to push
 * @param kmer_vars 	constant variables needed for kmer construction
 */

void mm_push_kmer(void *km, mm128_v *p, int glbl_index, kmer_info *kmer_vars)
{
	assert(glbl_index >= 0 && glbl_index < kmer_vars->len); // safeguard

	uint64_t shift1 = 2 * (kmer_vars->k - 1);
	uint64_t mask = (1ULL<<2*kmer_vars->k) - 1; 
	uint64_t kmer[2] = {0,0};
	int kmer_span = 0;
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

		kv_push(mm128_t, km, *p, info);
	}
}