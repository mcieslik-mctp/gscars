#include <Rcpp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>
#include <zlib.h>
#include <stdint.h>
using namespace Rcpp;


extern "C" 
{
  
  #include "barcodes.h"
  #include "preprocess.h"
  #include "seqtk/kseq.h"
  KSEQ_INIT(gzFile, gzread)
  
}

#define BC_LEN		16
#define MATE1_TRIM	7

//' Count barcodes
//'
//' @param inp_fn input FQ file
//' @param out_fn output counts file
//' @param wl_fn output counts file
//' @export
// [[Rcpp::export]]
StringVector countBarcodes(std::string inp_fn, std::string out_fn, std::string wl_fn) {

  FILE *inp_file = fopen(inp_fn.c_str(), "r");
  if (inp_file == NULL) {
    stop("inp_fn could not be opened");
  }

  FILE *out_file = fopen(out_fn.c_str(), "wb");
  if (out_file == NULL) {
    stop("out_fn could not be opened");
  }

  FILE *wl_file = fopen(wl_fn.c_str(), "r");
  if (wl_file == NULL) {
    stop("wl_fn could not be opened");
  }
  
  BarcodeDict wldict;
  wl_read(&wldict, wl_file);
  count_barcodes(&wldict, inp_file);
  wl_serialize(&wldict, out_file);

  /* clean-up */
  fclose(inp_file);
  fclose(out_file);
  fclose(wl_file);
  wl_dealloc(&wldict);
  return out_fn;
}


//' Preprocess barcoded FASTQ files.
//'
//' @param cts_inp_fn input FQ file read 1
//' @param fq1_inp_fn input FQ file read 1
//' @param fq2_inp_fn input FQ file read 2
//' @param fq1_out_fn output FQ file read 1
//' @param fq2_out_fn output FQ file read 2
//' @export
// [[Rcpp::export]]
StringVector preprocessFastq(std::string cts_inp_fn, std::string fq1_inp_fn, std::string fq2_inp_fn,
                                                     std::string fq1_out_fn, std::string fq2_out_fn) {
  
  if (access(cts_inp_fn.c_str(), F_OK) == -1) {
    stop("cts_inp_fn could not be opened for reading");
  }
  
  if(
     (access(fq1_inp_fn.c_str(), F_OK) == -1) ||
     (access(fq2_inp_fn.c_str(), F_OK) == -1)
     ) {
    stop("fq1_inp_fn or fq2_inp_fn could not be opened for reading");
  }

  FILE *fq_out[2];
  fq_out[0] = fopen(fq1_out_fn.c_str(), "wb");
  fq_out[1] = fopen(fq2_out_fn.c_str(), "wb");

  if (fq_out[0] == NULL || fq_out[1] == NULL) {
    stop("fq1_out_fn or fq2_out_fn could not be opened for writing");
  }

  // read counts
  FILE *cts_file = fopen(cts_inp_fn.c_str(), "rb");
  BarcodeDict wl;
  wl_deserialize(&wl, cts_file);
  wl_compute_priors(&wl);
  
  // open file handles
  gzFile fq_inp[2];
  fq_inp[0] = gzopen(fq1_inp_fn.c_str(), "rb");
  fq_inp[1] = gzopen(fq2_inp_fn.c_str(), "rb");

  kseq_t *ks_inp[2];
  ks_inp[0] = kseq_init(fq_inp[0]);
  ks_inp[1] = kseq_init(fq_inp[1]);

  kseq_t *r1, *r2;
  FILE *f1, *f2;

  char barcode[BC_LEN+1] = {'\0'};
  char barcode_qual[BC_LEN+1] = {'\0'};
  while ((kseq_read(ks_inp[0]) >= 0) && (kseq_read(ks_inp[1]) >= 0)) {

    

    
    f1 = fq_out[0];
    r1 = ks_inp[0];
      
    memcpy(barcode, r1->seq.s, BC_LEN);
    memcpy(barcode_qual, r1->qual.s, BC_LEN);
      
    fputc('>', f1);
    fputs(r1->name.s, f1);
    fputs(" RX:Z:", f1);
    fputs(barcode, f1);

    const int good_barcode = correct_barcode(barcode, barcode_qual, &wl);
    fputs(" BX:Z:", f1);
    fputs(barcode, f1);

    fputc('\n', f1);

      
    // fputs(barcode, f);
    // fputc('\n', f);
    // fputs("+", f);
    // fputc('\n', f);
    // fputs(barcode_qual, f);
    // fputc('\n', f);
      
  }
  
  // cleanup
  fclose(cts_file);  
  kseq_destroy(ks_inp[0]);
  kseq_destroy(ks_inp[1]);
  gzclose(fq_inp[0]);
  gzclose(fq_inp[1]);
  fclose(fq_out[0]);
  fclose(fq_out[1]);

  // output
  StringVector res(2);
  res[0] = fq1_out_fn;
  res[1] = fq2_out_fn;
  return res;
}


// void count_barcodes(BarcodeDict *bcdict, FILE *fq)
// {

// 	char barcode[BC_LEN+1];
// 	char id1[BUF_SIZE];
// 	char read1[BUF_SIZE];
// 	char sep1[BUF_SIZE];
// 	char qual1[BUF_SIZE];
// 	barcode[BC_LEN] = '\0';

//         char *ret;
// 	while (fgets(id1, BUF_SIZE, fq)) {
//                ret = fgets(read1, BUF_SIZE, fq);
//                assert(ret != NULL);
// 	       ret = fgets(sep1, BUF_SIZE, fq);
//                assert(ret != NULL);
// 	       ret = fgets(qual1, BUF_SIZE, fq);
//                assert(ret != NULL);

// 		for (size_t i = 0; i < BC_LEN; i++) {
// 			if (IS_ACGT(read1[i])) {
// 				barcode[i] = read1[i];
// 			} else {
// 				barcode[0] = '\0';
// 				break;
// 			}
// 		}
// 		if (barcode[0] != '\0') {
// 			const bc_t bc = encode_bc(barcode);
// 			wl_increment(bcdict, bc);
// 		}
// 	}
// }


// void preprocess_fastqs(const char *cts, const char *inp, const char *out)

// {

//   BarcodeDict wl;
//   assert(cts != NULL);
//   FILE *cts_file = fopen(cts, "rb");
//   if (cts_file == NULL) IOERROR(cts);
//   wl_deserialize(&wl, cts_file);
//   fclose(cts_file);
//   wl_compute_priors(&wl);
  
  
//   FILE *inp_file = (inp==NULL || EQ(inp, "-")) ? stdin : fopen(inp, "r");
//   if (inp_file == NULL) {
//     IOERROR(inp);
//   }
  
//   FILE *out_file = (out==NULL) ? stdout : fopen(out, "w");
//   if (out_file == NULL) {
//     IOERROR(out);
//   }

//   char barcode[BC_LEN+1];
//   char barcode_qual[BC_LEN+1];
//   char id1[BUF_SIZE];
//   char read1[BUF_SIZE];
//   char sep1[BUF_SIZE];
//   char qual1[BUF_SIZE];

//   barcode[BC_LEN] = '\0';
//   barcode_qual[BC_LEN] = '\0';

//   char *ret;
//   while (fgets(id1, BUF_SIZE, inp_file)) {
//          ret = fgets(read1, BUF_SIZE, inp_file);
//          assert(ret != NULL);
//          ret = fgets(sep1, BUF_SIZE, inp_file);
//          assert(ret != NULL);
//          ret = fgets(qual1, BUF_SIZE, inp_file);
//          assert(ret != NULL);

//     trim_after_space(id1);

//     size_t id1_len = strlen(id1);
//     size_t read1_len = strlen(read1);
//     size_t qual1_len = strlen(qual1);

//     assert(read1_len == qual1_len);
//     assert(read1_len > (BC_LEN + MATE1_TRIM));

    
//     memcpy(barcode, read1, BC_LEN);
//     memcpy(barcode_qual, qual1, BC_LEN);

//     /* trim off barcode and first MATE1_TRIM bases */
//     memmove(read1, &read1[BC_LEN + MATE1_TRIM], read1_len - (BC_LEN + MATE1_TRIM) + 1);
//     memmove(qual1, &qual1[BC_LEN + MATE1_TRIM], qual1_len - (BC_LEN + MATE1_TRIM) + 1);

//     const int good_barcode = correct_barcode(barcode, barcode_qual, &wl);

//     fputs(id1, out_file);
//     fputs(barcode, out_file);
//     fputs("\n", out_file);
//     fputs(read1, out_file);
//     fputs(sep1, out_file);
//     fputs(qual1, out_file);
//   }

//   wl_dealloc(&wl);
//   fclose(inp_file);
//   fclose(out_file);
// }
