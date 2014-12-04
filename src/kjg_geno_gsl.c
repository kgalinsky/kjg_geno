#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "kjg_geno.h"
#include "kjg_geno_gsl.h"

size_t KJG_GENO_GSL_ROWS = 256;

size_t
kjg_geno_gsl_slice (const kjg_geno* X, const size_t i, gsl_matrix* Y)
{
  return (kjg_geno_get_rows_normalized (X, i, Y->size1, Y->data));
}

void
kjg_geno_gsl_XTXA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B,
                   gsl_matrix *C)
{
  size_t i, r;                                                // row index
  gsl_matrix* Y = gsl_matrix_alloc (KJG_GENO_GSL_ROWS, X->n);

  gsl_matrix_set_zero (C);

  for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS)
    {
      r = kjg_geno_gsl_slice (X, i, Y);

      gsl_matrix_const_view Yi = gsl_matrix_const_submatrix (Y, 0, 0, r,
                                                             Y->size2);
      gsl_matrix_view Bi = gsl_matrix_submatrix (B, i, 0, r, B->size2);

      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, &Yi.matrix, A, 0,
                      &Bi.matrix);
      gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1, &Yi.matrix, &Bi.matrix, 1,
                      C);
    }

  gsl_matrix_free (Y);
}

void
kjg_geno_gsl_XA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B)
{
  size_t i, r;
  gsl_matrix* Y = gsl_matrix_alloc (KJG_GENO_GSL_ROWS, X->n);

  gsl_matrix_set_zero (B);

  for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS)
    {
      r = kjg_geno_gsl_slice (X, i, Y);

      gsl_matrix_const_view Yi = gsl_matrix_const_submatrix (Y, 0, 0, r,
                                                             Y->size2);
      gsl_matrix_view Bi = gsl_matrix_submatrix (B, i, 0, r, B->size2);

      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, &Yi.matrix, A, 0,
                      &Bi.matrix);
    }

  gsl_matrix_free (Y);
}

void
kjg_geno_gsl_XTA (const kjg_geno *X, const gsl_matrix *A, gsl_matrix *B)
{
  size_t i, r;
  gsl_matrix* Y = gsl_matrix_alloc (KJG_GENO_GSL_ROWS, X->n);

  gsl_matrix_set_zero (B);

  for (i = 0; i < X->m; i += KJG_GENO_GSL_ROWS)
    {
      r = kjg_geno_gsl_slice (X, i, Y);
      gsl_matrix_const_view Yi = gsl_matrix_const_submatrix (Y, 0, 0, r,
                                                             Y->size2);
      gsl_matrix_const_view Ai = gsl_matrix_const_submatrix (A, i, 0, r,
                                                             A->size2);
      gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1, &Yi.matrix, &Ai.matrix, 1,
                      B);
    }

  gsl_matrix_free (Y);
}
