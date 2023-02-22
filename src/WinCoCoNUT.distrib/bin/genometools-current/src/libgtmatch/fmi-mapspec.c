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
#include <stdint.h>
#include <inttypes.h>
#include "libgtcore/unused.h"
#include "fmindex.h"
#include "mapspec-def.h"
#define PRIu32
#define PRIu64
static void assignfmmapspecification(ArrayMapspecification *mapspectable,
                                     void *voidinfo,
                                     UNUSED bool writemode)
{
  Fmindexwithoptions *fmwithoptions = (Fmindexwithoptions *) voidinfo;
  Fmindex *fmindex;
  Mapspecification *mapspecptr;

  fmindex = fmwithoptions->fmptr;
  NEWMAPSPEC(fmindex->tfreq,Seqpos,(unsigned long) TFREQSIZE(fmindex->mapsize));
  NEWMAPSPEC(fmindex->superbfreq,Seqpos,
             (unsigned long) SUPERBFREQSIZE(fmindex->mapsize,
                                            fmindex->nofsuperblocks));
  NEWMAPSPEC(fmindex->markpostable,Seqpos,
             fmwithoptions->storeindexpos
             ? (unsigned long) MARKPOSTABLELENGTH(fmindex->bwtlength,
                                                  fmindex->markdist)
             : 0);
  NEWMAPSPEC(fmindex->boundarray,Bwtbound,(unsigned long) fmindex->numofcodes);
  NEWMAPSPEC(fmindex->specpos.spacePairBwtidx,PairBwtidx,
             fmwithoptions->storeindexpos
             ? fmindex->specpos.nextfreePairBwtidx
             : 0);
  NEWMAPSPEC(fmindex->bfreq,Uchar,
             (unsigned long) BFREQSIZE(fmindex->mapsize,fmindex->nofblocks));
}

int flushfmindex2file(FILE *fp,
                      Fmindex *fmindex,
                      bool storeindexpos,
                      Error *err)
{
  Fmindexwithoptions fmwithoptions;

  error_check(err);
  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return flushtheindex2file(fp,assignfmmapspecification,
                            (void *) &fmwithoptions,fmindex->sizeofindex,err);
}

int fillfmmapspecstartptr(Fmindex *fmindex,
                          bool storeindexpos,
                          const Str *tmpfilename,
                          Error *err)
{
  Fmindexwithoptions fmwithoptions;

  error_check(err);
  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return fillmapspecstartptr(assignfmmapspecification,
                             &fmindex->mappedptr,
                             (void *) &fmwithoptions,
                             tmpfilename,
                             fmindex->sizeofindex,
                             err);
}
