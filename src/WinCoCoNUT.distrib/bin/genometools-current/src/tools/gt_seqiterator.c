/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include <sys/types.h>
#include <sys/stat.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"
#include "libgtcore/progressbar.h"
#include "libgtcore/disc_distri.h"
#include "libgtmatch/format64.h"
#include "libgtmatch/stamp.h"
#include "tools/gt_seqiterator.h"

#define BUCKETSIZE 100
#define PRIu64
typedef struct
{
  bool verbose,
       dodistlen,
       doastretch;
} Seqiteroptions;

static OPrval parse_options(Seqiteroptions *seqiteroptions,
                            int *parsed_args,int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *optionverbose, *optiondistlen, *optionastretch;
  OPrval oprval;

  error_check(err);

  op = option_parser_new("[options] file [...]",
                         "Parse the supplied Fasta files.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  optionverbose = option_new_bool("v","be verbose",
                                  &seqiteroptions->verbose,false);
  option_parser_add_option(op, optionverbose);

  optiondistlen = option_new_bool("distlen",
                                  "show distribution of sequence length",
                                  &seqiteroptions->dodistlen,false);
  option_parser_add_option(op, optiondistlen);

  optionastretch = option_new_bool("astretch",
                                   "show distribution of A-substrings",
                                   &seqiteroptions->doastretch,false);
  option_exclude(optiondistlen, optionastretch);
  option_parser_add_option(op, optionastretch);

  option_parser_set_min_args(op, 1U);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void showdistseqlen(unsigned long key, unsigned long long value,
                           UNUSED void *data)
{
  unsigned long distvalue;

  distvalue = (unsigned long) value; /* XXX: is this cast always OK? */
  printf("%lu--%lu %lu\n",
         BUCKETSIZE * key,
         BUCKETSIZE * (key+1) - 1,
         distvalue);
}

typedef struct
{
  unsigned long long sumA, *mmercount;
  unsigned long maxvalue, minkey;
} Astretchinfo;

static void showastretches(unsigned long key, unsigned long long value,
                           void *data)
{
  Astretchinfo *astretchinfo = (Astretchinfo *) data;

  astretchinfo->sumA += value * (unsigned long long) key;
  if (key > astretchinfo->maxvalue)
  {
    astretchinfo->maxvalue = key;
  }
  /*@ignore@*/
  printf("%lu %llu\n", key,value);
  /*@end@*/
}

static void showmeroccurrence(unsigned long key, unsigned long long value,
                              void *data)
{
  unsigned long len;
  Astretchinfo *astretchinfo = (Astretchinfo *) data;

  for (len=astretchinfo->minkey; len<= key; len++)
  {
    assert(len <= astretchinfo->maxvalue);
    astretchinfo->mmercount[len] += value * (unsigned long long) (key-len+1);
  }
}

static unsigned long long accumulateastretch(DiscDistri *distastretch,
                                             const Uchar *sequence,
                                             unsigned long len)
{
  unsigned long i, lenofastretch = 0;
  unsigned long long countA = 0;

  for (i=0; i<len; i++)
  {
    if (sequence[i] == 'A' || sequence[i] == 'a')
    {
      countA++;
      lenofastretch++;
    } else
    {
      if (lenofastretch > 0)
      {
        disc_distri_add(distastretch,lenofastretch);
        lenofastretch = 0;
      }
    }
  }
  if (lenofastretch > 0)
  {
    disc_distri_add(distastretch,lenofastretch);
  }
  return countA;
}

static void processastretches(const DiscDistri *distastretch,
                              UNUSED unsigned long long countA)
{
  Astretchinfo astretchinfo;
  unsigned long len;

  astretchinfo.sumA = 0;
  astretchinfo.maxvalue = 0;
  astretchinfo.minkey = 10UL;
  disc_distri_foreach(distastretch,showastretches,&astretchinfo);
  astretchinfo.mmercount = ma_malloc(sizeof(*astretchinfo.mmercount) *
                                    (astretchinfo.maxvalue+1));
  memset(astretchinfo.mmercount,0,sizeof (*astretchinfo.mmercount) *
                                  (astretchinfo.maxvalue+1));
  disc_distri_foreach(distastretch,showmeroccurrence,&astretchinfo);
  for (len=astretchinfo.minkey; len<=astretchinfo.maxvalue; len++)
  {
    /*@ignore@*/
    printf("a^{%lu} occurs %llu times\n", len,astretchinfo.mmercount[len]);
    /*@end@*/
  }
  assert(astretchinfo.sumA == countA);
  ma_free(astretchinfo.mmercount);
}

int gt_seqiterator(int argc, const char **argv, Error *err)
{
  StrArray *files;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  int i, parsed_args, had_err;
  off_t totalsize;
  DiscDistri *distseqlen = NULL;
  DiscDistri *distastretch = NULL;
  uint64_t numofseq = 0, sumlength = 0;
  unsigned long minlength = 0, maxlength = 0;
  unsigned long long countA = 0;
  bool minlengthdefined = false;
  Seqiteroptions seqiteroptions;

  error_check(err);

  /* option parsing */
  switch (parse_options(&seqiteroptions,&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
        return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
        return 0;
  }

  files = strarray_new();
  for (i = parsed_args; i < argc; i++)
  {
    strarray_add_cstr(files, argv[i]);
  }
  totalsize = files_estimate_total_size(files);
  printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  seqit = seqiterator_new(files, NULL, true);
  if (seqiteroptions.dodistlen)
  {
    distseqlen = disc_distri_new();
  }
  if (seqiteroptions.doastretch)
  {
    distastretch = disc_distri_new();
  }
  if (seqiteroptions.verbose)
  {
    progressbar_start(seqiterator_getcurrentcounter(seqit, (unsigned long long)
                                                           totalsize),
                                                           (unsigned long long)
                                                           totalsize);
  }
  while (true)
  {
    desc = NULL;
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, err);
    if (seqiteroptions.dodistlen)
    {
      if (!minlengthdefined || minlength > len)
      {
        minlength = len;
        minlengthdefined = true;
      }
      if (maxlength < len)
      {
        maxlength = len;
      }
      sumlength += (uint64_t) len;
      numofseq++;
      disc_distri_add(distseqlen,len/BUCKETSIZE);
    }
    if (seqiteroptions.doastretch)
    {
      countA += accumulateastretch(distastretch,sequence,len);
    }
    ma_free(desc);
    if (had_err != 1)
    {
      break;
    }
  }
  if (seqiteroptions.verbose)
  {
    progressbar_stop();
  }
  seqiterator_delete(seqit);
  strarray_delete(files);
  if (seqiteroptions.dodistlen)
  {
    printf("# " Formatuint64_t " sequences of average length %.2f\n",
             PRINTuint64_tcast(numofseq),(double) sumlength/numofseq);
    printf("# minimum length %lu\n",minlength);
    printf("# maximum length %lu\n",maxlength);
    printf("# distribution of sequence length in buckets of size %d\n",
           BUCKETSIZE);
    disc_distri_foreach(distseqlen,showdistseqlen,NULL);
    disc_distri_delete(distseqlen);
  }
  if (seqiteroptions.doastretch)
  {
    processastretches(distastretch,countA);
    disc_distri_delete(distastretch);
  }
  return had_err;
}
