/*
 * kjg_geno_rand.h
 *
 *  Created on: Nov 23, 2014
 *      Author: Kevin Galinsky
 */

#ifndef KJG_GENO_RAND_H_
#define KJG_GENO_RAND_H_

#include <gsl/gsl_rng.h>

#include "kjg_geno.h"

/**
 * kjg_geno_rand_star - simulate random star-shaped population
 * @param r gsl_rng
 * @param M SNPs
 * @param N samples per population
 * @param P populations
 * @param FST
 * @param MAF minor allele frequency range (0-0.5)
 * @return genotype matrix
 */

kjg_geno*
kjg_geno_rand_star (gsl_rng* r, const size_t M, const size_t N, const size_t P,
                    const double FST, const double *MAF);

/**
 * kjg_geno_rand_anc - random ancestral AF within MAF range
 * @param r gsl_rng
 * @param MAF[2] minor allele frequency range (0-0.5)
 * @return
 */

double
kjg_geno_rand_anc (gsl_rng* r, const double* MAF);

/**
 * kjg_geno_rand_star_AF - allele frequencies for star-shaped ancestry
 * @param r gsl_rng
 * @param AF[P] allele frequencies to generate
 * @param anc ancestral AF
 * @param FST
 * @param P populations
 */

void
kjg_geno_rand_star_AF (gsl_rng* r, double* AF, const double anc,
                       const double FST, const size_t P);

/**
 * kjg_geno_rand_row - simulate a random row
 * @param r gsl_rng
 * @param x[N*P] row to simulate
 * @param N individuals per subpopulation
 * @param P subpopulations
 * @param AF[P] allele frequencies
 */

void
kjg_geno_rand_row (gsl_rng* r, uint8_t* x, const size_t N, const size_t P,
                   const double* AF);

/**
 * kjg_geno_rand_ld_row - simulate a random row in LD with another
 * @param r gsl_rng
 * @param x[N*P] row to simulate
 * @param N individuals per subpopulation
 * @param P subpopulations
 * @param AF[P] allele frequencies
 * @param y[N*P] original genotype row
 * @param cor[P] correlations
 */

void
kjg_geno_rand_ld_row (gsl_rng* r, uint8_t* x, const size_t N, const size_t P,
                      const double* AF, const uint8_t *y, const double* cor);

#endif /* KJG_GENO_RAND_H_ */
