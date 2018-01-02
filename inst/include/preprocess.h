#ifndef PREPROCESS_H
#define PREPROCESS_H

/* corrects barcodes and generate new FASTQs */
/* This is essentially a translation of Long Ranger's barcode correction scheme. */
void preprocess_fastqs(const char *cts, const char *inp, const char *out);

/* performs initial barcode count */
void count_barcodes(BarcodeDict *bcdict, FILE *fq);

#endif /* PREPROCESS_H */
