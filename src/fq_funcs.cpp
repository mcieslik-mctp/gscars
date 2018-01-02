#include <Rcpp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>
using namespace Rcpp;


extern "C" 
{

 #include "barcodes.h"
 #include "preprocess.h"

}



//' Extract putative DNA molecules from barcoded BAM file.
//'
//' @param inp_file input FQ file
//' @param out_file output counts file
//' @param out_file output counts file
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
