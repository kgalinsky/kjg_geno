#ifndef KJG_GENO_GSL_H_
#define KJG_GENO_GSL_H_

#include <gsl/gsl_matrix.h>

#include "kjg_geno.h"

/**
 * Multiplies B=X*A1 and A2 = XT*B = XT*X*A1
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A1 some matrix
 * @param *B intermediate matrix
 * @param *A2 next matrix
 */

void
kjg_geno_gsl_XTXA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B,
                   gsl_matrix *C);

/**
 * Multiplies B = X*A
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A some matrix
 * @param *B another matrix
 */

void
kjg_geno_gsl_XA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B);

/**
 * Multiplies B = XT*A
 * @param X compressed genotype matrix
 * @param *M array of SNP means
 * @param *A some matrix
 * @param *B another matrix
 */

void
kjg_geno_gsl_XTA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B);

#endif /* KJG_GENO_GSL_H_ */
