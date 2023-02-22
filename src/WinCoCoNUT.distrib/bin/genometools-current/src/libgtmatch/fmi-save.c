/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include "libgtcore/fa.h"
#include "fmindex.h"

#include "fmi-mapspec.pr"
#include "opensfxfile.pr"
#define PRIu32
static int writefmascii (const Str *indexname,
                         const Fmindex *fm,
                         bool storeindexpos,
                         Error *err)
{
  FILE *fmafp;

  error_check(err);
  if ((fmafp = opensfxfile (indexname, FMASCIIFILESUFFIX,"wb",err)) == NULL)
  {
    return -1;
  }
  fprintf (fmafp, "bwtlength=" FormatSeqpos "\n",
           PRINTSeqposcast(fm->bwtlength));
  fprintf (fmafp, "longest=" FormatSeqpos "\n",
                   PRINTSeqposcast(fm->longestsuffixpos));
  fprintf (fmafp, "storeindexpos=%d\n", storeindexpos ? 1 : 0);
  fprintf (fmafp, "log2blocksize=%u\n", fm->log2bsize);
  fprintf (fmafp, "log2markdist=%u\n", fm->log2markdist);
  fprintf (fmafp, "specialcharacters=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.specialcharacters));
  fprintf (fmafp, "specialranges=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.specialranges));
  fprintf (fmafp, "realspecialranges=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.realspecialranges));
  fprintf (fmafp, "lengthofspecialprefix=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.lengthofspecialprefix));
  fprintf (fmafp, "lengthofspecialsuffix=" FormatSeqpos "\n",
                  PRINTSeqposcast(fm->specialcharinfo.lengthofspecialsuffix));
  fprintf (fmafp, "suffixlength=%u\n", fm->suffixlength);
  fa_xfclose(fmafp);
  return 0;
}

static int writefmdata (const Str *indexname,
                        Fmindex *fm,
                        bool storeindexpos,
                        Error *err)
{
  FILE *fp;

  error_check(err);
  if ((fp = opensfxfile (indexname, FMDATAFILESUFFIX,"wb",err)) == NULL)
  {
    return -1;
  }
  if (flushfmindex2file(fp,fm,storeindexpos,err) != 0)
  {
    return -2;
  }
  fa_xfclose(fp);
  return 0;
}

int saveFmindex (const Str *indexname,Fmindex *fm,
                 bool storeindexpos,Error *err)
{
  error_check(err);
  if (writefmascii (indexname, fm, storeindexpos,err) != 0)
  {
    return -1;
  }
  if (writefmdata (indexname, fm, storeindexpos,err) != 0)
  {
    return -2;
  }
  return 0;
}
