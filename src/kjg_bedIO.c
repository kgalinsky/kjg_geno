#include "kjg_bedIO.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kjg_util.h"

kjg_bedIO*
kjg_bedIO_fopen (const char* path, const char* mode, const size_t m,
                 const size_t n)
{
  if (mode[0] != 'r')
    return (NULL); // TODO support writing

  FILE* stream = fopen (path, mode);
  if (stream == NULL)
    return (NULL);

  kjg_bedIO pre =
    { m, n, ((n - 1) / 4 + 1), stream };
  kjg_bedIO* bp = malloc (sizeof(kjg_bedIO));
  memcpy (bp, &pre, sizeof(kjg_bedIO));

  char magic[3];
  fread (magic, 1, 3, stream);
  if ((magic[0] != 0x6c) || (magic[1] != 0x1b) || (magic[2] != 0x01))
    {
      fprintf (stderr, "Bad magic numbers: %02x %02x %02x\n", magic[0],
               magic[1], magic[2]);
      exit (1);
    }
  return (bp);
}

kjg_bedIO*
kjg_bedIO_bfile_fopen (const char* path, const char* mode)
{
  size_t m = 0, n = 0;

  char c;

    {
      FILE* stream = kjg_util_fopen_suffix (path, "bim", "r");
      if (stream == NULL)
        return (NULL);
      while ((c = fgetc (stream)) != EOF)
        if (c == '\n')
          m++;
      fclose (stream);
    }

    {
      FILE* stream = kjg_util_fopen_suffix (path, "bim", "r");

      if (stream == NULL)
        return (NULL);
      while ((c = fgetc (stream)) != EOF)
        if (c == '\n')
          n++;
      fclose (stream);
    }

  char* filename;
  asprintf (&filename, "%s.bed", path);
  kjg_bedIO* bp = kjg_bedIO_fopen (filename, mode, m, n);
  free (filename);

  return (bp);
}

int
kjg_bedIO_fclose (kjg_bedIO* bp)
{
  int r = fclose (bp->stream);
  free (bp);
  return (r);
}

kjg_geno*
kjg_bedIO_fread_geno (kjg_bedIO* bp, const uint8_t* SNPmask,
                      const uint8_t* indmask)
{
  if ((!SNPmask) && (!indmask))
    {
      kjg_geno* g = kjg_geno_alloc (bp->m, bp->n);
      fread (g->data, 1, g->tda * g->m, bp->stream);
      return (g);
    }

  size_t m = bp->m, n = bp->n;
  size_t i;

  if (SNPmask)
    for (i = 0; i < bp->m; i++)
      m -= SNPmask[i];

  if (indmask)
    for (i = 0; i < bp->n; i++)
      n -= indmask[i];

  kjg_geno* g = kjg_geno_alloc (m, n);

  uint8_t *data = g->data;
  if (indmask)
    {
      uint8_t *buffer = malloc (bp->tda);
      for (i = 0; i < m; i++)
        {
          if (SNPmask && SNPmask[i])
            fseek (bp->stream, bp->tda, SEEK_CUR);
          fread (buffer, 1, bp->tda, bp->stream);
          kjg_geno_repack (bp->n, indmask, buffer, data);
          data += g->tda;
        }
      free (buffer);
    }
  else
    {
      for (i = 0; i < m; i++)
        {
          if (SNPmask[i])
            fseek (bp->stream, bp->tda, SEEK_CUR);
          fread (data, 1, bp->tda, bp->stream);
          data += g->tda;

        }
    }
  return (g);
}
