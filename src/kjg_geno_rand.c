/*
 * kjg_geno_rand.c
 *
 *  Created on: Nov 23, 2014
 *      Author: Kevin Galinsky
 */

#include "kjg_geno_rand.h"

#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kjg_geno.h"

kjg_geno*
kjg_geno_rand_star (gsl_rng* r, const size_t M, const size_t N, const size_t P,
                    const double FST, const double* MAF)
{
  kjg_geno *X = kjg_geno_alloc (M, N * P);
  double* AF = malloc (P * sizeof(double));
  uint8_t *x = malloc (N * P * sizeof(uint8_t));

  size_t m;
  double F = (1 - FST) / FST;

  for (m = 0; m < M; m++)
    {
      double anc = kjg_geno_rand_anc (r, MAF);
      kjg_geno_rand_star_AF (r, AF, anc, F, P);
      kjg_geno_rand_row (r, x, N, P, AF);
      kjg_geno_set_row (X, m, x);
    }

  free (x);
  free (AF);
  return (X);
}

double
kjg_geno_rand_anc (gsl_rng* r, const double* MAF)
{
  double anc = gsl_ran_flat (r, MAF[0], MAF[1]);
  if (gsl_ran_bernoulli (r, 0.5))
    anc = 1 - anc;
  return (anc);
}

void
kjg_geno_rand_star_AF (gsl_rng* r, double* AF, const double anc, const double F,
                       const size_t P)
{
  double a = anc * F;
  double b = (1 - anc) * F;
  size_t p;
  for (p = 0; p < P; p++)
    AF[p] = gsl_ran_beta (r, a, b);
}

void
kjg_geno_rand_row (gsl_rng* r, uint8_t* x, const size_t N, const size_t P,
                   const double* AF)
{
  size_t p, n, i = 0;
  for (p = 0; p < P; p++)
    {
      double cut1 = AF[p] * AF[p];
      double cut2 = 1 - (1 - AF[p]) * (1 - AF[p]);
      for (n = 0; n < N; n++)
        {
          double u = gsl_rng_uniform (r);
          x[i++] = u < cut1 ? 0 : u < cut2 ? 2 : 3;
        }
    }
}

void
kjg_geno_rand_ld_row (gsl_rng* r, uint8_t* x, const size_t N, const size_t P,
                      const double* AF, const uint8_t *y, const double* cor)
{
  size_t p, n, i = 0;
  for (p = 0; p < P; p++)
    {
      double cut1 = AF[p] * AF[p];
      double cut2 = 1 - (1 - AF[p]) * (1 - AF[p]);
      for (n = 0; n < N; n++)
        {
          double u = gsl_rng_uniform (r);
          if (u < cor[p])
            x[i] = y[i++];
          else
            {
              u = gsl_rng_uniform (r);
              x[i++] = u < cut1 ? 0 : u < cut2 ? 2 : 3;
            }
        }
    }
}

