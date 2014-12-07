#ifndef KJG_BEDIO_H_
#define KJG_BEDIO_H_

#include <stdio.h>

#include "kjg_geno.h"

typedef struct
{
  const size_t m;   // number of SNPs
  const size_t n;   // number of samples
  const size_t tda;
  FILE* stream;
} kjg_bedIO;

/**
 * Opens a bed file
 * @param path path to file
 * @param mode mode to open
 * @param m number of SNPs
 * @param n number of individuals
 * @return pointer to kjg_bedIO struct
 */

kjg_bedIO*
kjg_bedIO_fopen (const char* path, const char* mode, const size_t m,
                 const size_t n);

kjg_bedIO*
kjg_bedIO_bfile_fopen (const char* path, const char* mode);

/**
 * Closes a bed file
 * @param *bp pointer to kjg_bedIO struct
 */

int
kjg_bedIO_fclose (kjg_bedIO* bp);

kjg_geno*
kjg_bedIO_fread_geno (kjg_bedIO* bp, const uint8_t* SNPmask,
                      const uint8_t* indmask);

#endif /* KJG_BEDIO_H_ */
