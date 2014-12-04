/*
 * kjg_geno.c
 *
 *  Created on: Jul 31, 2013
 *      Author: kjg063
 *
 * Functions for working with genotype data.
 */

#include "kjg_geno.h"

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

// Set up lookup variables

#define P1(a, b, c) \
    { KJG_GENO_PACK(a, b, c, 0), \
      KJG_GENO_PACK(a, b, c, 1), \
      KJG_GENO_PACK(a, b, c, 2), \
      KJG_GENO_PACK(a, b, c, 3) }
#define P2(a, b) { P1(a, b, 0), P1(a, b, 1), P1(a, b, 2), P1(a, b, 3) }
#define P3(a)    { P2(a, 0),    P2(a, 1),    P2(a, 2),    P2(a, 3)    }
const uint8_t KJG_GENO_PACK_LOOKUP[4][4][4][4] =
  { P3(0), P3(1), P3(2), P3(3) };

#define U1(p) \
    KJG_GENO_UNPACK(p), \
    KJG_GENO_UNPACK(p + 1), \
    KJG_GENO_UNPACK(p + 2), \
    KJG_GENO_UNPACK(p + 3)
#define U2(p) U1(p), U1(p +  4), U1(p +  8), U1(p + 12)
#define U3(p) U2(p), U2(p + 16), U2(p + 32), U2(p + 48)
const uint8_t KJG_GENO_UNPACK_LOOKUP[256][4] =
  { U3(0), U3(64), U3(128),
U3(192) };

#define A1(p) \
    KJG_GENO_SUM_ALT(p), \
    KJG_GENO_SUM_ALT(p + 1), \
    KJG_GENO_SUM_ALT(p + 2), \
    KJG_GENO_SUM_ALT(p + 3)
#define A2(p) A1(p), A1(p +  4), A1(p +  8), A1(p + 12)
#define A3(p) A2(p), A2(p + 16), A2(p + 32), A2(p + 48)
const uint8_t KJG_GENO_SUM_ALT_LOOKUP[256] =
  { A3(0), A3(64), A3(128), A3(192) };

#define C1(p) \
    KJG_GENO_COUNT(p), \
    KJG_GENO_COUNT(p + 1), \
    KJG_GENO_COUNT(p + 2), \
    KJG_GENO_COUNT(p + 3)
#define C2(p) C1(p), C1(p +  4), C1(p +  8), C1(p + 12)
#define C3(p) C2(p), C2(p + 16), C2(p + 32), C2(p + 48)
const uint8_t KJG_GENO_COUNT_LOOKUP[256] =
  { C3(0), C3(64), C3(128), C3(192) };

// Functional interface

void
kjg_geno_pack (const size_t n, const uint8_t* u, uint8_t* p)
{
  size_t i = 0, j = 0;

  // pack the whole chunks
  for (; i < n - 4; (i += 4), j++)
    p[j] = kjg_geno_pack_unit (&u[i]);

  // pack the last chunk
  uint8_t remainder[4] =
    { 0, 0, 0, 0 };
  for (; i < n; i++)
    remainder[i % 4] = u[i];

  p[j] = kjg_geno_pack_unit (remainder);
}

void
kjg_geno_unpack (const size_t n, const uint8_t* p, uint8_t* u)
{
  size_t i = 0, j = 0;
  for (; i < n - 4; (i += 4), j++)
    memcpy (&u[i], KJG_GENO_UNPACK_LOOKUP[p[j]], 4);

  memcpy (&u[i], KJG_GENO_UNPACK_LOOKUP[p[j]], n - i);
}

size_t
kjg_geno_repack (const size_t n, const uint8_t* mask, const uint8_t* p1,
                 uint8_t* p2)
{
  size_t i, j, k; // iterate through p1
  size_t a = 0, b = 0; // iterate through p2
  uint8_t u2[4];

  for (i = 0, j = 0; i < n - 4; (i += 4), j++)
    {
      if (mask[i] && mask[i + 1] && mask[i + 2] && mask[i + 3])
        continue;

      const uint8_t *u1 = KJG_GENO_UNPACK_LOOKUP[p1[j]];
      for (k = 0; k < 4; k++)
        {
          if (!mask[i + k])
            {
              u2[b++] = u1[k];

              if (b == 4)
                {
                  p2[a++] = kjg_geno_pack_unit (u2);
                  b = 0;
                }
            }
        }
    }

  const uint8_t *u1 = KJG_GENO_UNPACK_LOOKUP[p1[j]];
  for (; i < n; i++)
    {
      if (!mask[i])
        {
          u2[b++] = u1[i % 4];

          if (b == 4)
            {
              p2[a++] = kjg_geno_pack_unit (u2);
              b = 0;
            }
        }
    }

  if (b)
    {
      for (; b < 4; b++)
        {
          u2[b] = 0;
        }
      p2[a] = kjg_geno_pack_unit (u2);
    }

  return (a);
}

size_t
kjg_geno_repack4 (const size_t n, const uint8_t* mask, const uint8_t* p1,
                  uint8_t* p2)
{
  size_t i = 0, j = 0, k = 0;
  for (; i < n; (i += 4), j++)
    if (!mask[j])
      p2[k++] = p1[j];

  return (k);
}

size_t
kjg_geno_sum_alt (const size_t n, const uint8_t* p)
{
  size_t i = 0, j = 0, a = 0;

  for (; i < n; (i += 4), j++)
    a += KJG_GENO_SUM_ALT_LOOKUP[p[j]];

  return (a);
}

size_t
kjg_geno_count (const size_t n, const uint8_t* p)
{
  size_t i = 0, j = 0, c = 0;

  for (; i < n; (i += 4), j++)
    c += KJG_GENO_COUNT_LOOKUP[p[j]];

  c -= (n + 4 - i) % 4;

  return (c);
}

double
kjg_geno_af (const size_t n, const uint8_t* p)
{
  return (((double) kjg_geno_sum_alt (n, p)) / kjg_geno_count (n, p) / 2);
}

// Normalization

int
kjg_geno_norm (const double p, double s[4])
{
  if (p == 0 || p == 1)
    {             // check for homogeneous population
      memset (s, 0, 4);                // zero out scaling array
      return (1);                     // error
    }

  size_t i;

  double d = sqrt (2 * p * (1 - p));   // Var(G) = 2pq
  s[0] = -2 * p / d;
  s[2] = (1 - 2 * p) / d;
  s[3] = (2 - 2 * p) / d;
  s[1] = 0;

  return (0);
}

// Constructor/Destructor

kjg_geno*
kjg_geno_alloc (size_t m, size_t n)
{
  kjg_geno pre =
    { m, n, (((n - 1) / 4) + 1) };
  kjg_geno* g = malloc (sizeof(kjg_geno));

  memcpy (g, &pre, sizeof(kjg_geno));

  g->data = malloc (sizeof(uint8_t) * m * g->tda);
  g->af = 0;
  g->norm = 0;

  return (g);
}

void
kjg_geno_free (kjg_geno* g)
{
  free (g->data);
  if (g->af != 0)
    free (g->af);
  if (g->norm != 0)
    free (g->norm);
  free (g);
}

// Getter/Setter

void
kjg_geno_get_row (const kjg_geno* g, const size_t i, uint8_t* x)
{
  kjg_geno_unpack (g->n, g->data + g->tda * i, x);
}

void
kjg_geno_get_row_normalized (const kjg_geno* g, const size_t i, double* y)
{
  const uint8_t* p = g->data + i * g->tda;
  const double* s = g->norm + i * 4;

  if (s[0] == 0)
    {
      memset (y, 0, sizeof(double) * g->n);
      return;
    }

  size_t j = 0, k = 0;
  for (; j < g->n - 4; (j += 4), k++)
    {
      const uint8_t* u = KJG_GENO_UNPACK_LOOKUP[p[k]];
      y[j] = s[u[0]];
      y[j + 1] = s[u[1]];
      y[j + 2] = s[u[2]];
      y[j + 3] = s[u[3]];
    }

  const uint8_t* u = KJG_GENO_UNPACK_LOOKUP[p[k]];
  for (; j < g->n; j++)
    y[j] = s[u[j % 4]];
}

size_t
kjg_geno_get_rows_normalized (const kjg_geno* g, const size_t i, const size_t r,
                              double* Y)
{
  size_t j;
  for (j = i; (j < i + r) && (j < g->m); j++, Y += g->n)
    kjg_geno_get_row_normalized (g, j, Y);
  return (j - i);
}

void
kjg_geno_set_row (kjg_geno* g, const size_t i, const uint8_t* x)
{
  kjg_geno_pack (g->n, x, g->data + g->tda * i);
}

void
kjg_geno_set_af (kjg_geno* g, double* af)
{
  if (af != 0)
    {
      g->af = af;
      return;
    }

  if (g->af == 0)
    g->af = malloc (sizeof(double) * g->m);

  size_t i;
  for (i = 0; i < g->m; i++)
    g->af[i] = kjg_geno_af (g->n, g->data + (i * g->tda));
}

void
kjg_geno_set_norm (kjg_geno* g, double* norm)
{
  if (norm != 0)
    {
      g->norm = norm;
      return;
    }

  if (g->af == 0)
    kjg_geno_set_af (g, 0);
  if (g->norm == 0)
    g->norm = malloc (sizeof(double) * g->m * 4);

  size_t i;
  for (i = 0; i < g->m; i++)
    kjg_geno_norm (g->af[i], g->norm + i * 4);
}
