#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "barcodes.h"
#include "preprocess.h"

/* lengths */
#define BC_LEN		16
#define MATE1_TRIM	7
#define BUF_SIZE	1024

void count_barcodes(BarcodeDict *bcdict, FILE *fq)
{

	char barcode[BC_LEN+1];
	char id1[BUF_SIZE];
	char read1[BUF_SIZE];
	char sep1[BUF_SIZE];
	char qual1[BUF_SIZE];
	barcode[BC_LEN] = '\0';

        char *ret;
	while (fgets(id1, BUF_SIZE, fq)) {
               ret = fgets(read1, BUF_SIZE, fq);
               assert(ret != NULL);
	       ret = fgets(sep1, BUF_SIZE, fq);
               assert(ret != NULL);
	       ret = fgets(qual1, BUF_SIZE, fq);
               assert(ret != NULL);

		for (size_t i = 0; i < BC_LEN; i++) {
			if (IS_ACGT(read1[i])) {
				barcode[i] = read1[i];
			} else {
				barcode[0] = '\0';
				break;
			}
		}

		if (barcode[0] != '\0') {
			const bc_t bc = encode_bc(barcode);
			wl_increment(bcdict, bc);
		}
	}

}

static int correct_barcode(char *barcode, char *barcode_qual, BarcodeDict *wl)
{
#define ILLUMINA_QUAL_OFFSET 33

	uint8_t quals[BC_LEN];

	int n_count = 0;
	for (size_t i = 0; i < BC_LEN; i++) {
		barcode[i] = toupper(barcode[i]);
		quals[i] = barcode_qual[i] - ILLUMINA_QUAL_OFFSET;

		if (!IS_ACGT(barcode[i])) {
			++n_count;
		}
	}

	const bc_t bc0 = ((n_count == 0) ? encode_bc(barcode) : 0);
	BarcodeInfo *bcinfo = ((n_count == 0) ? wl_lookup(wl, bc0) : NULL);

#define CAND_BUF_SIZE (120 * 4 * 4)
	struct bc_string { char bc_str[BC_LEN]; };
	struct bc_string bc_cands[CAND_BUF_SIZE];
	double bc_cand_probs[CAND_BUF_SIZE];
	size_t n_cands = 0;
#undef CAND_BUF_SIZE

	if (bcinfo == NULL) {
		if (n_count > 1) {
			return 0;
		}

		/* examine Hamming-1 neighbors */
		for (size_t i = 0; i < BC_LEN; i++) {
			const char prev = barcode[i];

			if (n_count > 0 && IS_ACGT(prev))
				continue;

			for (size_t j = 0; j < 4; j++) {
				const char new = "ACGT"[j];

				if (new == prev)
					continue;

				barcode[i] = new;
				const bc_t bc = encode_bc(barcode);
				bcinfo = wl_lookup(wl, bc);

				if (bcinfo != NULL) {
					const double prior = bcinfo->prior;
					const double edit_log_prob = MIN(33.0, (double)(quals[i]));
					const double edit_prob = pow(10.0, (-edit_log_prob/10.0));
					const double p = prior*edit_prob;

					memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
					bc_cand_probs[n_cands] = p;
					++n_cands;
				}
			}

			barcode[i] = prev;
		}
	} else {
		assert(n_count == 0);
		memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
		bc_cand_probs[n_cands] = bcinfo->prior;
		++n_cands;

		/* examine Hamming-2 neighbors */
		for (size_t i1 = 0; i1 < BC_LEN; i1++) {
			const char prev1 = barcode[i1];

			for (size_t j1 = 0; j1 < 4; j1++) {
				const char new1 = "ACGT"[j1];

				if (new1 == prev1)
					continue;

				barcode[i1] = new1;

				for (size_t i2 = i1 + 1; i2 < BC_LEN; i2++) {
					const char prev2 = barcode[i2];

					for (size_t j2 = 0; j2 < 4; j2++) {
						const char new2 = "ACGT"[j2];

						if (new2 == prev2)
							continue;

						barcode[i2] = new2;
						const bc_t bc = encode_bc(barcode);
						bcinfo = wl_lookup(wl, bc);

						if (bcinfo != NULL) {
							const double prior = bcinfo->prior;

							const double e1 = MAX(3.0, (double)(quals[i1]) - 1.0);
							const double e2 = MAX(3.0, (double)(quals[i2]) - 1.0);

							const double edit1_log_prob = MIN(33.0, e1);
							const double edit2_log_prob = MIN(33.0, e2);

							const double edit_prob = pow(10.0, (-edit1_log_prob/10.0)) * pow(10.0, (-edit2_log_prob/10.0));
							const double p = prior*edit_prob;

							memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
							bc_cand_probs[n_cands] = p;
							++n_cands;
						}
					}

					barcode[i2] = prev2;
				}

				barcode[i1] = prev1;
			}
		}
	}

	if (n_cands > 0) {
		double total = bc_cand_probs[0];
		size_t max = 0;

		for (size_t i = 1; i < n_cands; i++) {
			total += bc_cand_probs[i];

			if (bc_cand_probs[i] > bc_cand_probs[max])
				max = i;
		}

		if (bc_cand_probs[max]/total > BC_CONF_THRESH) {
			memcpy(barcode, bc_cands[max].bc_str, BC_LEN);
			return 1;
		}
	}

	return 0;

#undef ILLUMINA_QUAL_OFFSET
}

void preprocess_fastqs(const char *cts, const char *inp, const char *out)

{

  BarcodeDict wl;
  assert(cts != NULL);
  FILE *cts_file = fopen(cts, "rb");
  if (cts_file == NULL) IOERROR(cts);
  wl_deserialize(&wl, cts_file);
  fclose(cts_file);
  wl_compute_priors(&wl);
  
  
  FILE *inp_file = (inp==NULL || EQ(inp, "-")) ? stdin : fopen(inp, "r");
  if (inp_file == NULL) {
    IOERROR(inp);
  }
  
  FILE *out_file = (out==NULL) ? stdout : fopen(out, "w");
  if (out_file == NULL) {
    IOERROR(out);
  }

  char barcode[BC_LEN+1];
  char barcode_qual[BC_LEN+1];
  char id1[BUF_SIZE];
  char read1[BUF_SIZE];
  char sep1[BUF_SIZE];
  char qual1[BUF_SIZE];

  barcode[BC_LEN] = '\0';
  barcode_qual[BC_LEN] = '\0';

  char *ret;
  while (fgets(id1, BUF_SIZE, inp_file)) {
         ret = fgets(read1, BUF_SIZE, inp_file);
         assert(ret != NULL);
         ret = fgets(sep1, BUF_SIZE, inp_file);
         assert(ret != NULL);
         ret = fgets(qual1, BUF_SIZE, inp_file);
         assert(ret != NULL);

    trim_after_space(id1);

    size_t id1_len = strlen(id1);
    size_t read1_len = strlen(read1);
    size_t qual1_len = strlen(qual1);

    assert(read1_len == qual1_len);
    assert(read1_len > (BC_LEN + MATE1_TRIM));

    
    memcpy(barcode, read1, BC_LEN);
    memcpy(barcode_qual, qual1, BC_LEN);

    /* trim off barcode and first MATE1_TRIM bases */
    memmove(read1, &read1[BC_LEN + MATE1_TRIM], read1_len - (BC_LEN + MATE1_TRIM) + 1);
    memmove(qual1, &qual1[BC_LEN + MATE1_TRIM], qual1_len - (BC_LEN + MATE1_TRIM) + 1);

    const int good_barcode = correct_barcode(barcode, barcode_qual, &wl);

    fputs(id1, out_file);
    fputs(barcode, out_file);
    fputs("\n", out_file);
    fputs(read1, out_file);
    fputs(sep1, out_file);
    fputs(qual1, out_file);
  }

  wl_dealloc(&wl);
  fclose(inp_file);
  fclose(out_file);
}
