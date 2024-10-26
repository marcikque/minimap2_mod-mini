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


typedef struct {		// struct to simplify mm_push_kmer's arguments (improves readability and reduces duplicate code)
	void *km;			// thread-local memory pool
	mm128_v *p;

	mm128_t *min;		// pointer to current minimal element
	int *min_pos;		// pointer to variable holding pos of current minimizer within tmer buffer

	mm128_t *buf;		// ring buffer containing tmers
	int *buf_pos;		// pointer to variable holding current pos in tmer's ring buffer

	mm128_t *k_buf;		// ring buffer containing kmers
	int *k_buf_pos;		// pointer to variable holding current pos in kmer's ring buffer

	int *prev_push;		// position of last pushed kmer within the sequence
	int len;			// length of DNA sequence
} kmer_info;

int static mm_push_kmer(int w, int k, int t, kmer_info *kmer_vars);

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
	int r = (int) ceil(log2(w+k-1) / 2) + 1; //default: w = 10, k = 15, r = 4
	int t = r + ((k-r) % w);

	uint64_t shift1 = 2 * (t - 1);
	uint64_t mask = (1ULL<<2*t) - 1;
	uint64_t tmer[2] = {0,0};
	int l = 0; // bases processed since last reset/start
	int buf_pos = 0;
	int min_pos = 0;
	int tmer_span = 0; // bases included in the current tmer (default t, changes under HPC)
	int prev_push = -1; // saves last pushed kmer's position within the sequence
	mm128_t buf[256];
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	// variables needed for simultaneous kmer creation
	uint64_t k_shift1 = 2 * (k - 1);
	uint64_t k_mask = (1ULL<<2*k) - 1;
	uint64_t kmer[2] = {0,0};
	int k_buf_pos = 0;
	int kmer_span = 0;
	mm128_t k_buf[256];
	tiny_queue_t k_tq;

	// creating the necessary structs to call mm_push_kmer
	kmer_info kmer_vars = {
		.km = km,
		.p = p,
		.min = &min,
		.min_pos = &min_pos,
		.buf = buf,
		.buf_pos = &buf_pos,	
		.k_buf = k_buf,
		.k_buf_pos = &k_buf_pos,
		.prev_push = &prev_push,
		.len = len
	};

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56bits for kmers, 8bits for metadata

	memset(buf, 0xff, (w + k - t) * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	memset(k_buf, 0xff, w * 16);
	memset(&k_tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (int i = 0; i < len; ++i) {
		mm128_t info = { UINT64_MAX, UINT64_MAX }; 
		mm128_t k_info = { UINT64_MAX, UINT64_MAX }; 
		int c = seq_nt4_table[(uint8_t)str[i]];

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
				tq_push(&k_tq, skip_len);

				tmer_span += skip_len; // count number of bases within tmer, can be >t
				kmer_span += skip_len;
				if (tq.count > t) {
					tmer_span -= tq_shift(&tq);
				} 
				if (k_tq.count > k) {
					kmer_span -= tq_shift(&k_tq);
				} 
			} else {
				tmer_span = l + 1 < t? l + 1 : t;
				kmer_span = l + 1 < k? l + 1 : k;
			}

			tmer[0] = (tmer[0] << 2 | c) & mask;           // forward t-mer: remove first base, append new base
			tmer[1] = (tmer[1] >> 2) | (3ULL^c) << shift1; // reverse t-mer: remove last base, prepend complement of new base (A<->T, C<->G)

			kmer[0] = (kmer[0] << 2 | c) & k_mask;           // forward k-mer: remove first base, append new base
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << k_shift1; // reverse k-mer: remove last base, prepend complement of new base (A<->T, C<->G)


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


			if (kmer[0] == kmer[1]) {
				continue; // symmetric tmer, strand unknown
			}
			int k_z = kmer[0] < kmer[1]? 0 : 1; // strand: forward or backward (lexicographically smaller is chosen)

			if (l >= k && kmer_span < 256) {
				k_info.x = hash64(kmer[k_z], k_mask) << 8 | kmer_span; 
				// kmer.x[63:8] = kmer[k_z] hashed, trimmed
				// kmer.x[7:0]  = kmer_span
				k_info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | k_z;
				// kmer.y[63:32] = sequence ref. ID
				// kmer.y[31:1]  = last position
				// kmer.y[0]     = strand
			}
		
		// resetting, essentially splitting the sequence, at unambigious bases
		} else {
			l = 0;
			tq.count = 0;
			tq.front = 0;
			tmer_span = 0;
			buf_pos = 0;

			k_tq.count = 0;
			k_tq.front = 0;
			kmer_span = 0;
			k_buf_pos = 0;

			min_pos = 0;
			min = (mm128_t){ UINT64_MAX, UINT64_MAX };
		}

		buf[buf_pos] = info;
		k_buf[k_buf_pos] = k_info;

		// case: new minimum found
		if (info.x < min.x) {
			min = info;
			min_pos = buf_pos;
			
			if (l >= w + k) { // check if we fully processed the first window
				prev_push = mm_push_kmer(w, k, t, &kmer_vars);
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

			prev_push = mm_push_kmer(w, k, t, &kmer_vars);
			
		// case: consecutive windows share minimizing tmer
		} else {
			if (l >= w + k) { // check if we fully processed the first window
				prev_push = mm_push_kmer(w, k, t, &kmer_vars);
			}
		}

		// case: 1st window is fully processed, not partial anymore, we need to push its minimizer
		if (l == w + k - 1) {
			prev_push = mm_push_kmer(w, k, t, &kmer_vars);
		}

		//inc tmer buffer position
		if (++buf_pos == w + k - t) {
			buf_pos = 0;
		}

		//inc kmer buffer position
		if (++k_buf_pos == w) {
			k_buf_pos = 0;
		}
	}
}


/** 
 * Finds the kmer corresponding to the current minimal tmer in mod-minimizer fashion.
 * If it has not been pushed yet, we push it (duplicate protection).
 * Returns the kmer's last position for the next iteration's duplication check.
 * @param w 		window length
 * @param k 		kmer length
 * @param t 		tmer length
 * @param kmer_vars 	struct to simplify the arguments
*/
int static mm_push_kmer(int w, int k, int t, kmer_info *kmer_vars) {
	int buf_size = w + k - t;
    int k_buf_size = w;

	// 1. calculate the minimal element's positon within the window
	int window_pos = (*kmer_vars->min_pos - *kmer_vars->buf_pos + (buf_size - 1)) % buf_size;
	// 2. mod-minimizer: modulo operation
	window_pos %= w;
	// 3. calculate the kmer buffer index of the window position
	int k_buf_idx = (*kmer_vars->k_buf_pos + 1 + window_pos) % k_buf_size;
	// 4. retrieve the corresponding kmer
	mm128_t k_target = kmer_vars->k_buf[k_buf_idx];

    // retrieve the target's end position within the sequence
    uint32_t k_target_glbl_pos = k_target.y;
    k_target_glbl_pos >>= 1;
    assert(k_target_glbl_pos >= 0 && k_target_glbl_pos < kmer_vars->len);

    // Push the k-mer if not already pushed
    if (*kmer_vars->prev_push != k_target_glbl_pos) {
        kv_push(mm128_t, kmer_vars->km, *kmer_vars->p, k_target);
    }

    return k_target_glbl_pos;
}