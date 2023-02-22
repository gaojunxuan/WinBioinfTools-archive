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
#include "libgtcore/chardef.h"
#include "libgtcore/fa.h"
#include "emimergeesa.h"
#include "esafileend.h"
#include "fmindex.h"
#include "divmodmul.h"
#include "spacedef.h"

#include "esa-merge.pr"
#include "esa-map.pr"
#include "encseq2offset.pr"
#include "opensfxfile.pr"
#include "mkidxcpy.pr"
#include "fmi-keyval.pr"
#include "fmi-mapspec.pr"
#define PRIu32
 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Seqpos);

static void allocatefmtables(Fmindex *fm,bool storeindexpos)
{
  ALLOCASSIGNSPACE (fm->tfreq, NULL, Seqpos,TFREQSIZE(fm->mapsize));
  ALLOCASSIGNSPACE (fm->superbfreq, NULL, Seqpos ,
                    SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks));
  if (storeindexpos)
  {
    ALLOCASSIGNSPACE (fm->markpostable,NULL,Seqpos,
                      MARKPOSTABLELENGTH(fm->bwtlength,fm->markdist));
    fm->specpos.nextfreePairBwtidx = 0;
    fm->specpos.allocatedPairBwtidx
      = (unsigned long) determinenumberofspecialstostore(&fm->specialcharinfo);
    printf("# %lu wildcards in the last " FormatSeqpos
           " characters (%.2f)\n",
            (unsigned long) (fm->specialcharinfo.specialcharacters -
                             fm->specpos.allocatedPairBwtidx),
            PRINTSeqposcast(fm->specialcharinfo.specialcharacters),
            (double) (fm->specialcharinfo.specialcharacters -
                      fm->specpos.allocatedPairBwtidx)/
                          fm->specialcharinfo.specialcharacters);
    ALLOCASSIGNSPACE(fm->specpos.spacePairBwtidx,NULL,PairBwtidx,
                     fm->specpos.allocatedPairBwtidx);
  } else
  {
    INITARRAY(&fm->specpos,PairBwtidx);
    fm->markpostable = NULL;
  }
  ALLOCASSIGNSPACE (fm->bfreq, NULL, Uchar,
                    BFREQSIZE(fm->mapsize,fm->nofblocks));
}

static void set0frequencies(Fmindex *fm)
{
  Seqpos i;

  for (i = 0; i < (Seqpos) TFREQSIZE(fm->mapsize); i++)
  {
    fm->tfreq[i] = 0;
  }
  for (i = 0; i < (Seqpos) BFREQSIZE(fm->mapsize,fm->nofblocks); i++)
  {
    fm->bfreq[i] = 0;
  }
  for (i = 0; i < (Seqpos) SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks); i++)
  {
    fm->superbfreq[i] = 0;
  }
}

static void finalizefmfrequencies(Fmindex *fm)
{
  unsigned int j;
  Seqpos i, *freqptr;

  for (j = 2U; j <= fm->mapsize; j++)
  {
    fm->tfreq[j] += fm->tfreq[j - 1];
  }
  freqptr = fm->superbfreq;
  for (j = 0; j < fm->mapsize; j++)
  {
    for (i = (Seqpos) 2; i < fm->nofsuperblocks; i++)
    {
      freqptr[i] += freqptr[i-1];
    }
    freqptr += fm->nofsuperblocks;
  }
}

static void showconstructionmessage(const Str *indexname,
                                    Seqpos totallength,
                                    unsigned long fmsize,
                                    unsigned int log2bsize,
                                    unsigned int log2markdist,
                                    unsigned int mapsize)
{
  printf("# construct fmindex \"%s\" for bsize=%lu, superbsize=%lu,",
          str_get(indexname),
          (unsigned long) POW2(log2bsize),
          (unsigned long) POW2(log2markdist));
  printf(" len=" FormatSeqpos ", alphasize=%u: size ",
          PRINTSeqposcast(totallength),
          mapsize-1);
  printf("%lu bytes, space overhead %.2f\n",
          fmsize,
          (double) fmsize/(double) (totallength+1));
}

static int nextesamergedsufbwttabvalues(DefinedSeqpos *longest,
                                        Uchar *bwtvalue,
                                        Seqpos *suftabvalue,
                                        Emissionmergedesa *emmesa,
                                        const Seqpos *sequenceoffsettable,
                                        Seqpos bwtpos,
                                        Error *err)
{
  Indexedsuffix indexedsuffix;

  error_check(err);
  if (emmesa->buf.nextaccessidx >= emmesa->buf.nextstoreidx)
  {
    if (emmesa->numofentries == 0)
    {
      return 0;
    }
    if (stepdeleteandinsertothersuffixes(emmesa,err) != 0)
    {
      return -1;
    }
    if (emmesa->buf.nextstoreidx == 0)
    {
      return 0;
    }
    emmesa->buf.nextaccessidx = 0;
  }
  indexedsuffix = emmesa->buf.suftabstore[emmesa->buf.nextaccessidx];
  *suftabvalue = sequenceoffsettable[indexedsuffix.idx] +
                 indexedsuffix.startpos;
  if (indexedsuffix.startpos == 0)
  {
    if (indexedsuffix.idx == 0)
    {
      if (longest->defined)
      {
        error_set(err,"longest is already defined as " FormatSeqpos,
                      longest->valueseqpos);
        return -2;
      }
      longest->defined = true;
      longest->valueseqpos = bwtpos;
      *bwtvalue = (Uchar) UNDEFBWTCHAR;
    } else
    {
      *bwtvalue = (Uchar) SEPARATOR;
    }
  } else
  {
    *bwtvalue
      = getencodedchar( /* Random access */
           emmesa->suffixarraytable[indexedsuffix.idx].encseq,
           indexedsuffix.startpos-1,
           emmesa->suffixarraytable[indexedsuffix.idx].readmode);
  }
  emmesa->buf.nextaccessidx++;
  return 1;
}

int sufbwt2fmindex(Fmindex *fmindex,
                   unsigned int log2bsize,
                   unsigned int log2markdist,
                   const Str *outfmindex,
                   const StrArray *indexnametab,
                   bool storeindexpos,
                   Verboseinfo *verboseinfo,
                   Error *err)
{
  Suffixarray suffixarray;
  Emissionmergedesa emmesa;
  Uchar cc;
  Seqpos bwtpos,
         totallength,
         suftabvalue = 0,
         *sequenceoffsettable = NULL,
         firstignorespecial = 0,
         nextmark,
         *markptr,
         nextprogress,
         tmpsuftabvalue,
         stepprogress;
  unsigned int mapsize = 0,
               suffixlength = 0,
               numofindexes;
  int retval;
  DefinedSeqpos longest = { false, 0 };
  PairBwtidx *pairptr;
  FILE *outbwt = NULL;
  Str *tmpfilename = NULL;
  Specialcharinfo specialcharinfo;
  bool haserr = false;

  error_check(err);
  longest.defined = false;
  longest.valueseqpos = 0;
  numofindexes = (unsigned int) strarray_size(indexnametab);
  if (numofindexes == 1U)
  {
    Str *indexname = strarray_get_str(indexnametab,0);

    if (streamsuffixarray(&suffixarray,
                         &totallength,
                         SARR_BWTTAB | (storeindexpos ? SARR_SUFTAB : 0),
                         indexname,
                         verboseinfo,
                         err) != 0)
    {
      haserr = true;
    }
    if (!haserr)
    {
      mapsize = getmapsizeAlphabet(suffixarray.alpha);
      specialcharinfo = suffixarray.specialcharinfo;
      firstignorespecial = totallength - specialcharinfo.specialcharacters;
      if (makeindexfilecopy(outfmindex,indexname,ALPHABETFILESUFFIX,0,err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (makeindexfilecopy(outfmindex,
                            indexname,
                            BWTTABSUFFIX,
                            firstignorespecial,
                            err) != 0)
      {
        haserr = true;
      }
    }
  } else
  {
    if (initEmissionmergedesa(&emmesa,
                             indexnametab,
                             SARR_ESQTAB | SARR_SUFTAB | SARR_LCPTAB,
                             verboseinfo,
                             err) != 0)
    {
      haserr = true;
    }
    if (!haserr)
    {
      Str *indexname = strarray_get_str(indexnametab,0);
      suffixlength = 0;
      if (makeindexfilecopy(outfmindex,indexname,ALPHABETFILESUFFIX,0,err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      sequenceoffsettable = encseqtable2seqoffsets(&totallength,
                                                   &specialcharinfo,
                                                   emmesa.suffixarraytable,
                                                   numofindexes);
      if (sequenceoffsettable == NULL)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      longest.defined = false;
      longest.valueseqpos = 0;
      outbwt = opensfxfile(outfmindex,BWTTABSUFFIX,"wb",err);
      if (outbwt == NULL)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      mapsize = getmapsizeAlphabet(emmesa.alpha);
      firstignorespecial = totallength - specialcharinfo.specialcharacters;
    }
  }
  if (!haserr)
  {
    printf("# firstignorespecial=" FormatSeqpos "\n",
              PRINTSeqposcast(firstignorespecial));
    computefmkeyvalues (fmindex,
                        totallength+1,
                        log2bsize,
                        log2markdist,
                        mapsize,
                        suffixlength,
                        storeindexpos,
                        &specialcharinfo);
    showconstructionmessage(outfmindex,
                            totallength,
                            fmindex->sizeofindex,
                            log2bsize,
                            log2markdist,
                            mapsize);
    allocatefmtables(fmindex,storeindexpos);
    set0frequencies(fmindex);
    if (storeindexpos)
    {
      markptr = fmindex->markpostable;
    } else
    {
      markptr = NULL;
    }
    nextprogress = stepprogress = totallength/78;
    for (bwtpos = 0, nextmark = 0; ; bwtpos++)
    {
      if (numofindexes == 1U)
      {
        if (storeindexpos)
        {
          retval = readnextSeqposfromstream(&tmpsuftabvalue,
                                            &suffixarray.suftabstream,
                                            err);
          if (retval < 0)
          {
            haserr = true;
            break;
          }
          if (retval == 0)
          {
            break;
          }
          suftabvalue = (Seqpos) tmpsuftabvalue;
        }
        retval = readnextUcharfromstream(&cc,&suffixarray.bwttabstream,err);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          break;
        }
      } else
      {
        retval = nextesamergedsufbwttabvalues(&longest,
                                              &cc,
                                              &suftabvalue,
                                              &emmesa,
                                              sequenceoffsettable,
                                              bwtpos,
                                              err);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          break;
        }
        if (fwrite(&cc,
                  sizeof (Uchar),
                  (size_t) 1,
                  outbwt) != (size_t) 1)
        {
          haserr = true;
          break;
        }
      }
      if (bwtpos == nextprogress)
      {
        if (bwtpos == stepprogress)
        {
          (void) putchar('#');
        }
        (void) putchar('.');
        (void) fflush(stdout);
        nextprogress += stepprogress;
      }
      if (storeindexpos && bwtpos == nextmark)
      {
        *markptr++ = suftabvalue;
        nextmark += fmindex->markdist;
      }
      if (ISBWTSPECIAL(cc))
      {
        if (storeindexpos && bwtpos < firstignorespecial)
        {
          pairptr = fmindex->specpos.spacePairBwtidx +
                    fmindex->specpos.nextfreePairBwtidx++;
          if (pairptr >= fmindex->specpos.spacePairBwtidx +
                         fmindex->specpos.allocatedPairBwtidx)
          {
            error_set(err,"program error: not enough space for specpos");
            haserr = true;
            break;
          }
          pairptr->bwtpos = bwtpos;
          pairptr->suftabvalue = suftabvalue;
        }
      } else
      {
        fmindex->tfreq[cc+1]++;
        fmindex->bfreq[(cc * fmindex->nofblocks) +
                       (bwtpos >> fmindex->log2bsize)]++;
        fmindex->superbfreq[(cc * fmindex->nofsuperblocks) +
                       (bwtpos >> fmindex->log2superbsize) + 1]++;
      }
    }
  }
  if (!haserr)
  {
    if (storeindexpos &&
        fmindex->specpos.allocatedPairBwtidx !=
        fmindex->specpos.nextfreePairBwtidx)
    {
      error_set(err,"program error: too much space for specpos: "
                    "allocated = %lu != %lu = used",
                    fmindex->specpos.allocatedPairBwtidx,
                    fmindex->specpos.nextfreePairBwtidx);
      haserr = true;
    }
  }
  if (!haserr)
  {
    (void) putchar('\n');
    finalizefmfrequencies(fmindex);
    if (fmindex->suffixlength > 0)
    {
      ALLOCASSIGNSPACE(fmindex->boundarray,NULL,Bwtbound,fmindex->numofcodes);
    }
    if (numofindexes == 1U)
    {
      fmindex->longestsuffixpos = suffixarray.longest.valueseqpos;
      freesuffixarray(&suffixarray);
    } else
    {
      if (!longest.defined)
      {
        error_set(err,"longest is not defined after merging");
        haserr = true;
      }
      if (!haserr)
      {
        fmindex->longestsuffixpos = longest.valueseqpos;
      }
      fa_xfclose(outbwt);
      wraptEmissionmergedesa(&emmesa);
    }
  }
  FREESPACE(sequenceoffsettable);
  FREESPACE(tmpfilename);
  return haserr ? -1 : 0;
}
