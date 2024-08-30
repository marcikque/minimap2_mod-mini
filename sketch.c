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

typedef struct {		// kmer selection information
	int k;				// kmer length
	int w;				// window size
	int t;				// tmer length
	int buf_pos;		// current ring-buffer position
	int w_relative;		// mod-sampled element's position relative to current window
	mm128_t min;		// min.x, min.y
} idx_info;

typedef struct {		// kmer construction information
	const char *str;	// sequence
	int len;			// sequence length
	uint32_t rid;		// reference ID
	int is_hpc;			// enable homopolymer-compression
} kmer_info;

void mm_push_kmer(void *km, mm128_v *p, mm128_t *buf, idx_info *idx_vars, kmer_info *kmer_vars);
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
	int r = (int) ceil(log2(w+k-1) / 2) + 1; // r=4 default | TODO: Threshold experiment
	int t = r + ((k-r) % w); // Section 4.2 (11:13), Def. 16
	assert(k > r); 

	uint64_t shift1 = 2 * (t - 1);
	uint64_t mask = (1ULL<<2*t) - 1;
	uint64_t tmer[2] = {0,0};
	int l = 0;
	int buf_pos = 0;
	int min_pos = 0;
	int tmer_span = 0;
	mm128_t buf[256];
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56bits for kmers, 8bits for metadata

	memset(buf, 0xff, (w + k - t) * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (int i = 0; i < len; ++i) {
		mm128_t info = { UINT64_MAX, UINT64_MAX }; 
			// -> saveguard to check for info.x,.y==MAX
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
		}

		buf[buf_pos] = info;

		/*
		- windows go from 0 to w-1, contain w kmers, the last of which starts at w-1 and ends at w-1+k-1
		- windows contain w+k-t tmers, the last of which starts at w+k-t-1 and ends at w+k-2
		- this code is constantly operating on a window's last tmer
		- buf_pos is our relative position within the ring buffer, it has no say about our relative position within the window
			(which is always at the very end of the very last tmer: position w-1+k-1)
		*/

		/*
		issue: due to modulo, some of the resulting indices point to buffer positions that contain UINT64_MAX
		result: invalid kmers and sequence accesses
		solution: instead of starting index of the kmer, we now transmit the tmer to check for invalid tmers
		*/

		// special case: 1st window identical kmers
		if (l == w + k - 1 && min.x != UINT64_MAX) {
			int w_relative = 0;
			for (int j = buf_pos + 1; j < w + k - t; ++j, ++w_relative) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					idx_info idx_vars = {
						.k = k,
						.w = w,
						.t = t,
						.buf_pos = buf_pos,
						.w_relative = w_relative,
						.min = min
					};
					kmer_info kmer_vars = {
						.str = str,
						.len = len,
						.rid = rid,
						.is_hpc = is_hpc
					};
					mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
				}
			}
			for (int j = 0; j < buf_pos; ++j, ++w_relative) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					idx_info idx_vars = {
						.k = k,
						.w = w,
						.t = t,
						.buf_pos = buf_pos,
						.w_relative = w_relative,
						.min = min
					};
					kmer_info kmer_vars = {
						.str = str,
						.len = len,
						.rid = rid,
						.is_hpc = is_hpc
					};
					mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
				}
			}	
			
		}

		// case: new minimum found, push the old one
		if (info.x <= min.x) { 
			if (l >= w + k && min.x != UINT64_MAX) {
				// the old min was the min of the previous window
				// operate relative to previous window -> decrement buf_pos and shift into valid index
				int prev_buf_pos = (buf_pos - 1 + w + k - t) % (w + k - t);
				int w_relative;
				if (prev_buf_pos >= min_pos) {
					w_relative = (w + k - t - 1 - prev_buf_pos + min_pos); 
				} else {
					w_relative = (min_pos - prev_buf_pos - 1);
				}
				idx_info idx_vars = {
					.k = k,
					.w = w,
					.t = t,
					.buf_pos = prev_buf_pos,
					.w_relative = w_relative,
					.min = min
				};
				kmer_info kmer_vars = {
					.str = str,
					.len = len,
					.rid = rid,
					.is_hpc = is_hpc
				};
				mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
			}
			min = info;
			min_pos = buf_pos;
			

		//case: old minimum moved out of the window (buffer)
		} else if (buf_pos == min_pos) {
			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				// old min had to be at index 0 of previous window
				idx_info idx_vars = {
					.k = k,
					.w = w,
					.t = t,
					.buf_pos = buf_pos,
					.w_relative = -1,
					.min = min
				};
				kmer_info kmer_vars = {
					.str = str,
					.len = len,
					.rid = rid,
					.is_hpc = is_hpc
				};
				mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
			}

			min.x = UINT64_MAX;
			for (int j = buf_pos + 1; j <  w + k - t; ++j) {
				if (min.x >= buf[j].x) {
					min = buf[j];
					min_pos = j;
				}
			}
			for (int j = 0; j <= buf_pos; ++j) {
				if (min.x >= buf[j].x) {
					min = buf[j];
					min_pos = j;
				}
			}

			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				// simply: the for loops ensure that we go through the window in order, starting at pos 0 (1 after buf_pos)
				int w_relative = 0;
				for (int j = buf_pos + 1; j < w + k - t; ++j, ++w_relative) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						idx_info idx_vars = {
							.k = k,
							.w = w,
							.t = t,
							.buf_pos = buf_pos,
							.w_relative = w_relative,
							.min = min
						};
						kmer_info kmer_vars = {
							.str = str,
							.len = len,
							.rid = rid,
							.is_hpc = is_hpc
						};
						mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
					}
				}
				for (int j = 0; j <= buf_pos; ++j, ++w_relative) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						idx_info idx_vars = {
							.k = k,
							.w = w,
							.t = t,
							.buf_pos = buf_pos,
							.w_relative = w_relative,
							.min = min
						};
						kmer_info kmer_vars = {
							.str = str,
							.len = len,
							.rid = rid,
							.is_hpc = is_hpc
						};
						mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
					}
				}
			}
		}
		if (++buf_pos == w + k - t) {
			buf_pos = 0;
		}
	}
	if (min.x != UINT64_MAX) { 
		// PRE: j>=m: X = (w + k - t - 1 - j + m) % w
		// POST: j<m: X = (m - j - 1) % w

		// Translate back: (X + j + 1) % (w+k-t)
		int w_relative;
		if (buf_pos >= min_pos) {
			w_relative = (w + k - t - 1 - buf_pos + min_pos); 
		} else {
			w_relative = (min_pos - buf_pos - 1);
		}
		idx_info idx_vars = {
			.k = k,
			.w = w,
			.t = t,
			.buf_pos = buf_pos,
			.w_relative = w_relative,
			.min = min
		};
		kmer_info kmer_vars = {
			.str = str,
			.len = len,
			.rid = rid,
			.is_hpc = is_hpc
		};
		mm_push_kmer(km, p, buf, &idx_vars, &kmer_vars);
	}
}


/**
 * Constructs and pushes the kmer starting at index i onto the list of minimizers p
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 * @param i		 index of kmer to push
 */

void mm_push_kmer(void *km, mm128_v *p, mm128_t *buf, idx_info *idx_vars, kmer_info *kmer_vars)
{
	mm128_t target_tmer;
	// special case for "old minimum moved out of the window (buffer)"
	if (idx_vars->w_relative == -1) {
		target_tmer = idx_vars->min;
	} else {
		int min_index = idx_vars->w_relative % idx_vars->w;
		int buf_relative = (idx_vars->buf_pos + 1 + min_index) % (idx_vars->w + idx_vars->k - idx_vars->t);
		target_tmer = buf[buf_relative];
	}

	// check needed as with mod-sampling we might try to push with invalid args
	if (target_tmer.y == UINT64_MAX || target_tmer.x == UINT64_MAX) {
		return;
	}
	int glbl_relative = target_tmer.y;

	assert(glbl_relative >= 0 && kmer_vars->len > glbl_relative);

	uint64_t shift1 = 2 * (idx_vars->k - 1);
	uint64_t mask = (1ULL<<2*idx_vars->k) - 1; 
	uint64_t kmer[2] = {0,0};
	int kmer_span = 0;
	int l = 0;
	int j;

	for (j=glbl_relative; j<kmer_vars->len && l<idx_vars->k; ++j) {
		int c = seq_nt4_table[(uint8_t)kmer_vars->str[j]];

		if (c>=4) { // ambiguous base
			return;
		}

		if (kmer_vars->is_hpc) { // TODO: dissect HPC functionality and potentially upgrade
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

	if (l >= idx_vars->k && kmer_span < 256) {
		

		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z = kmer[0] < kmer[1]? 0 : 1;

		info.x = hash64(kmer[z], mask) << 8 | kmer_span; 
		// kmer.x[63:8] = kmer[z] hashed, trimmed
		// kmer.x[7:0]  = kmer_span
		info.y = (uint64_t)kmer_vars->rid<<32 | (uint32_t)j<<1 | z; //FIXME: Logic inconsistency due to j (0) being erronerous
		
		//fprintf(stderr, "Checkstat: i = %d, j = %d, l = %d, len = %d\n", glbl_relative, j, l, kmer_vars->len);
		// kmer.y[63:31] = sequence ref. ID
		// kmer.y[7:1]   = last position
		// kmer.y[0]     = strand

		kv_push(mm128_t, km, *p, info);
	}
}