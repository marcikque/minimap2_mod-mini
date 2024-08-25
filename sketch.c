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

void mm_push_index(void *km, const char *str, int len, int k, uint32_t rid, int is_hpc, mm128_v *p, int i);

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
	assert(!is_hpc); // TODO: blocked for index translation testing
	
	// "The mod-minimizer: a simple and efficientsampling algorithm for long k-mers"
	int r = (int) ceil(log2(w+k-1) / 2) + 1; // r=4 default | TODO: Threshold experiment
	int t = r + ((k-r) % w); // Section 4.2 (11:13), Def. 16
	assert(k > r); 


	uint64_t shift1 = 2 * (t - 1);
	uint64_t mask = (1ULL<<2*t) - 1;
	uint64_t tmer[2] = {0,0};
	int i= 0;
	int j = 0;
	int l = 0;
	int buf_pos = 0;
	int min_pos = 0;
	int tmer_span = 0;
	mm128_t buf[256]; // ring buffer for tmers in current window
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56bits for kmers, 8bits for metadata

	memset(buf, 0xff, (w + k - t) * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {

		int c = seq_nt4_table[(uint8_t)str[i]];

		mm128_t info = { UINT64_MAX, UINT64_MAX };

		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				// next base same? -> loop till different
				// e.g. GTAAAAAAAACTG -> GTACTG
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

			tmer[0] = (tmer[0] << 2 | c) & mask;           // forward tmer
			tmer[1] = (tmer[1] >> 2) | (3ULL^c) << shift1; // reverse tmer: reverse & A<->T, C<->G
			z = tmer[0] < tmer[1]? 0 : 1; // strand: forward or backward (lexicographically smaller is chosen)
			++l;

			if (l >= t && tmer_span < 256) {
				info.x = hash64(tmer[z], mask) << 8 | tmer_span; 
				// tmer.x[63:8] = tmer[z] hashed, trimmed
				// tmer.x[7:0]  = tmer_span
				info.y = (uint64_t)(i - tmer_span + 1); 
				// tmer.y[63:0] = starting position
			}
		} else { // base no good -> invalid, reset and split sequence
			l = 0, tq.count = tq.front = 0, tmer_span = 0;
		}

		buf[buf_pos] = info;

		/*
		- windows go from 0 to w-1, contain w kmers, the last of which starts at w-1 and ends at w+k-2
		- windows contain w+k-t tmers, the last of which starts at w+k-t-1 and ends at w+k-2
		- this code is constantly operating on a window's last tmer
		*/
		

		// special case: 1st window identical tmers
		if (l == w + k - 1 && min.x != UINT64_MAX) { 
			for (j = buf_pos + 1; j < w + k - t; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					int min_index = i - (w + k - t) + (j % w);

					mm_push_index(km, str, len, k, rid, 0, p, min_index);
				}
			}
			/* 
			special case part 2, removed as an experiment since l==w+k-1 ==> buf_pos == 0 
			(except perhaps at c >= 4)
			for (j = 0; j < buf_pos; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					kv_push(mm128_t, km, *p, buf[j]);
				}
			} 
			*/
		}

		if (info.x <= min.x) {
			
			if (l >= w + k && min.x != UINT64_MAX) {
				int min_index = i - (w + k - t) + ((w + k - t - i) % w);
				mm_push_index(km, str, len, k, rid, 0, p, min_index);
			}
			min = info, min_pos = buf_pos;

		// case: old minimum moves out of the buffer (window)
		} else if (buf_pos == min_pos) {

			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				int min_index = i - (w + k - t - 2);
				mm_push_index(km, str, len, k, rid, 0, p, min_index);
			}

			// relative position = w + k - t - stepsBack
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) {
				if (min.x >= buf[j].x) {
					min = buf[j], min_pos = j;
				}
			}
			for (j = 0; j <= buf_pos; ++j){ 
				if (min.x >= buf[j].x) { 
					min = buf[j], min_pos = j;
				}
			}

			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				int w_rel = 0;
				for (j = buf_pos + 1; j < w + k - t; ++j, ++w_rel) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						int min_index = i - (w + k - t) + (w_rel % w);
						mm_push_index(km, str, len, k, rid, 0, p, min_index);
					}
				}
				for (j = 0; j <= buf_pos; ++j, ++w_rel) {
					if (min.x == buf[j].x && min.y != buf[j].y) {

						int min_index = i - (w + k - t) + (w_rel % w);
						mm_push_index(km, str, len, k, rid, 0, p, min_index);
					}
				}
			}
		}

		if (++buf_pos == w + k - t) {
			buf_pos = 0;
		}
	}

	if (min.x != UINT64_MAX) {
		int min_index = (i-(w + k - 1)) + ((w+k-1-(i-min.y))%w);
		mm_push_index(km, str, len, k, rid, 0, p, min_index);
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
void mm_push_index(void *km, const char *str, int len, int k, uint32_t rid, int is_hpc, mm128_v *p, int i) {
	assert(len > i);

	uint64_t shift1 = 2 * (k - 1);
	uint64_t mask = (1ULL<<2*k) - 1; 
	uint64_t kmer[2] = {0,0};
	int kmer_span = 0;
	int l = 0;
	int j;

	for (j=i; j<len && l<k; ++j) {
		int c = seq_nt4_table[(uint8_t)str[j]];

		if (c>=4) { // ambiguous base
			return;
		}

		if (is_hpc) { // TODO: dissect HPC functionality
			int skip_len = 1;
			if (j + 1 < len && seq_nt4_table[(uint8_t)str[j + 1]] == c) { 

				for (skip_len = 2; j + skip_len < len; ++skip_len) {
					if (seq_nt4_table[(uint8_t)str[j + skip_len]] != c) {
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

	if (l == k && kmer_span < 256) {
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z = kmer[0] < kmer[1]? 0 : 1;

		info.x = hash64(kmer[z], mask) << 8 | kmer_span; 
		// kmer.x[63:8] = kmer[z] hashed, trimmed
		// kmer.x[7:0]  = kmer_span

		info.y = (uint64_t)rid<<32 | (uint32_t)j<<1 | z;
		// kmer.y[63:31] = sequence ref. ID
		// kmer.y[7:1]   = last position
		// kmer.y[0]     = strand

		kv_push(mm128_t, km, *p, info);
	}
}