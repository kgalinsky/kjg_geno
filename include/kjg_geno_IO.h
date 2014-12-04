/*
 * @file kjg_geno_IO.h
 * @brief Reads geno files
 */

#ifndef KJG_GENO_IO_H_
#define KJG_GENO_IO_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "kjg_geno.h"

// Struct for reading geno files

typedef struct
{
  const size_t m;   // number of SNPs
  const size_t n;   // number of samples
  FILE* stream;
} kjg_geno_IO;

/**
 * Opens a geno file
 * @param *path path to file
 * @param *mode mode to open (supports only read for now)
 * @return pointer to kjg_geno_IO struct
 */

kjg_geno_IO*
kjg_geno_IO_fopen (const char* path, const char* mode);

/**
 * Closes a geno file
 * @param *gp pointer to kjg_geno_IO struct
 */

int
kjg_geno_IO_fclose (kjg_geno_IO* gp);

/**
 * Read geno file into struct
 *
 * @param *gp geno_IO struct
 * @return point to kjg_geno struct with data
 */

kjg_geno*
kjg_geno_IO_fread_geno (kjg_geno_IO* gp);

size_t
kjg_geno_IO_fread_chunk (kjg_geno_IO* gp, kjg_geno* g);

/**
 * Determine the number of individuals in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @return number of individuals
 */

size_t
kjg_geno_IO_num_ind (FILE* stream);

/**
 * Determine the number of SNPs in a geno file.
 *
 * @param *stream file pointer object for geno file
 * @param n number of individuals
 * @return number of SNPs
 */

size_t
kjg_geno_IO_num_snp (FILE* stream, size_t n);

/**
 * Convert a character buffer to geno.
 *
 * @param *buffer character buffer
 * @param *x genotype array
 * @param n number of individuals
 */

void
kjg_geno_IO_char2int (const char* buffer, uint8_t* x, const size_t n);

#endif /* KJG_GENO_IO_H_ */
