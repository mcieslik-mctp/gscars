#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "barcodes.h"
#include "preprocess.h"
#include "main.h"


static void print_help_and_exit(const char *argv0, int error)
{
#define P(...) fprintf(out, __VA_ARGS__)
  
	FILE *out = error ? stderr : stdout;
        P("fastq10x version %s\n", VERSION);
	P("usage: %s <count|preprocess|help>\n", argv0);
	exit(error ? EXIT_FAILURE : EXIT_SUCCESS);

#undef P
}

int main(const int argc, char *argv[])
{

  const char *argv0 = argv[0];

  if (argc < 2) {
    print_help_and_exit(argv0, 1);
  }
  
  const char *mode = argv[1];

  /* HELP */
  if (EQ(mode, "help")) {
    print_help_and_exit(argv0, 0);
  }

  /* COUNT */
  if (EQ(mode, "count")) {

    char c;
    char *wl = NULL;
    char *inp = NULL;
    char *out = NULL;
    while ((c = getopt(argc, argv, "w:i:o:")) != -1) {
      switch (c) {
      case 'w':
        wl = strdup(optarg);
        break;
      case 'i':
        inp = strdup(optarg);
        break;
      case 'o':
        out = strdup(optarg);
        break;
      default:
        print_help_and_exit(argv0, 1);
      }
    }

    if (wl == NULL) {
      fprintf(stderr, "error: specify barcode whitelist with -w\n");
      exit(EXIT_FAILURE);
    }

    FILE *inp_file = (inp==NULL || EQ(inp, "-")) ? stdin : fopen(inp, "r");
    if (inp_file == NULL) {
      IOERROR(inp);
    }

    FILE *out_file = (out==NULL) ? stdout : fopen(out, "wb");
    if (out_file == NULL) {
      IOERROR(out);
    }

    /* perform barcode counting */
    BarcodeDict wldict;
    wl_read(&wldict, wl);
    count_barcodes(&wldict, inp_file);
    wl_serialize(&wldict, out_file);

    /* clean-up */
    fclose(inp_file);
    fclose(out_file);
    wl_dealloc(&wldict);
    free(wl);
    free(inp);
    free(out);
    
    return EXIT_SUCCESS;
  }
  
  if (EQ(mode, "preproc")) {
          
    char c;
    char *cts = NULL;
    char *inp = NULL;
    char *out = NULL;

    while ((c = getopt(argc, argv, "c:i:o:")) != -1) {
      switch (c) {
      case 'c':
        cts = strdup(optarg);
        break;
      case 'i':
        inp = strdup(optarg);
        break;
      case 'o':
        out = strdup(optarg);
        break;
      default:
        print_help_and_exit(argv0, 1);
      }
    }
          
    if (cts == NULL) {
      fprintf(stderr, "error: specify barcode counts file with -c\n");
      exit(EXIT_FAILURE);
    }
          
    preprocess_fastqs(cts, inp, out);
    return EXIT_SUCCESS;
    
  }
        
  fprintf(stderr, "error: unrecognized mode\n");
  print_help_and_exit(argv0, 1);

}
