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

#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include "libgtcore/unused.h"
#include "sarr-def.h"
#include "spacedef.h"
#include "emimergeesa.h"
#include "encseq-def.h"
#include "merger-trie.h"
#include "verbose-def.h"

#include "esa-map.pr"
#define PRIu32
 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

static void fillandinsert(Mergertrierep *trierep,
                          unsigned int idx,
                          Seqpos suftabvalue,
                          Mergertrienode *node,
                          UNUSED uint64_t ident)
{
  Suffixinfo sinfo;

  sinfo.idx = idx;
  sinfo.startpos = suftabvalue;
#ifdef WITHTRIEIDENT
  sinfo.ident = ident;
#endif
  insertsuffixintomergertrie(trierep,node,&sinfo);
}

static int inputthesequences(Alphabet **alpha,
                             Seqpos *nextpostable,
                             Suffixarray *suffixarraytable,
                             const StrArray *indexnametab,
                             unsigned int demand,
                             Verboseinfo *verboseinfo,
                             Error *err)
{
  unsigned long idx;
  Str *indexname;
  Seqpos totallength;

  error_check(err);
  for (idx=0; idx<strarray_size(indexnametab); idx++)
  {
    indexname = strarray_get_str(indexnametab,idx);
    if (streamsuffixarray(&suffixarraytable[idx],
                          &totallength,
                          demand,
                          indexname,
                          verboseinfo,
                          err) != 0)
    {
      str_delete(indexname);
      return -1;
    }
    if (idx == 0)
    {
      *alpha = suffixarraytable[idx].alpha;
    }
    nextpostable[idx] = 0;
  }
  return 0;
}

static int insertfirstsuffixes(Mergertrierep *trierep,
                               Seqpos *nextpostable,
                               Suffixarray *suffixarraytable,
                               unsigned int numofindexes,
                               Error *err)
{
  unsigned int idx;
  Seqpos suftabvalue;
  int retval;

  error_check(err);
  for (idx=0; idx<numofindexes; idx++)
  {
    retval = readnextSeqposfromstream(&suftabvalue,
                                      &suffixarraytable[idx].suftabstream,
                                      err);
    if (retval < 0)
    {
      return -1;
    }
    if (retval == 0)
    {
      error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading suftab",__FILE__,__LINE__);
      return -2;
    }
    nextpostable[idx]++;
    fillandinsert(trierep,
                  idx,
                  suftabvalue,
                  trierep->root,
                  (uint64_t) idx);
  }
  return 0;
}

/*@null@*/ static Mergertrienode *findlargestnodeleqlcpvalue(
                                          Mergertrienode *smallest,
                                          Seqpos lcpvalue,
                                          Error *err)
{
  Mergertrienode *tmp;

  error_check(err);
  for (tmp = smallest->parent; tmp != NULL; tmp = tmp->parent)
  {
    if (tmp->depth <= lcpvalue)
    {
      return tmp;
    }
  }
  error_set(err,"path does not contain a node of depth <= " FormatSeqpos,
                PRINTSeqposcast(lcpvalue));
  return NULL;
}

int stepdeleteandinsertothersuffixes(Emissionmergedesa *emmesa, Error *err)
{
  Mergertrienode *tmpsmallestleaf, *tmplcpnode;
  Largelcpvalue tmpexception;
  Uchar tmpsmalllcpvalue;
  int retval;
  Seqpos tmpsuftabvalue,
         tmplcpvalue,
         tmplastbranchdepth;
  unsigned int tmpidx;

  error_check(err);
  for (emmesa->buf.nextstoreidx = 0;
      emmesa->numofentries > 0 &&
      emmesa->buf.nextstoreidx < (unsigned int) SIZEOFMERGERESULTBUFFER;
      emmesa->buf.nextstoreidx++)
  {
    tmpsmallestleaf = findsmallestnodeintrie(&emmesa->trierep);
    tmplastbranchdepth = tmpsmallestleaf->parent->depth;
    tmpidx = tmpsmallestleaf->suffixinfo.idx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].idx = tmpidx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].startpos
      = tmpsmallestleaf->suffixinfo.startpos;
    if (emmesa->nextpostable[tmpidx] >
       getencseqtotallength(emmesa->suffixarraytable[tmpidx].encseq))
    {
      deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
      emmesa->numofentries--;
    } else
    {
      retval = readnextUcharfromstream(&tmpsmalllcpvalue,
                                       &emmesa->suffixarraytable[tmpidx].
                                                lcptabstream,err);
      if (retval < 0)
      {
        return -1;
      }
      if (retval == 0)
      {
        error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading lcptab",__FILE__,__LINE__);
        return -2;
      }
      if (tmpsmalllcpvalue == (Uchar) UCHAR_MAX)
      {
        retval = readnextLargelcpvaluefromstream(
                               &tmpexception,
                               &emmesa->suffixarraytable[tmpidx].llvtabstream,
                               err);
        if (retval < 0)
        {
          return -3;
        }
        if (retval == 0)
        {
          error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading llvtab",__FILE__,__LINE__);
          return -4;
        }
        tmplcpvalue = tmpexception.value;
      } else
      {
        tmplcpvalue = (Seqpos) tmpsmalllcpvalue;
      }
      if (tmplcpvalue > tmplastbranchdepth)
      {
        tmplastbranchdepth = tmplcpvalue;
      }
      tmplcpnode = findlargestnodeleqlcpvalue(tmpsmallestleaf,tmplcpvalue,err);
      retval = readnextSeqposfromstream(&tmpsuftabvalue,
                                        &emmesa->suffixarraytable[tmpidx].
                                        suftabstream,
                                        err);
      if (retval < 0)
      {
        return -5;
      }
      if (retval == 0)
      {
        error_set(err,"file %s: line %d: unexpected end of file when "
                      "reading suftab",__FILE__,__LINE__);
        return -6;
      }
      emmesa->nextpostable[tmpidx]++;
      fillandinsert(&emmesa->trierep,
                    tmpidx,
                    tmpsuftabvalue,
                    tmplcpnode,
                    emmesa->ident++);
      tmpsmallestleaf = findsmallestnodeintrie(&emmesa->trierep);
      deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
    }
    if (emmesa->numofentries > 0)
    {
      emmesa->buf.lcptabstore[emmesa->buf.nextstoreidx]
        = tmplastbranchdepth;
      emmesa->buf.lastpage = false;
    } else
    {
      emmesa->buf.lastpage = true;
    }
  }
  return 0;
}

int initEmissionmergedesa(Emissionmergedesa *emmesa,
                          const StrArray *indexnametab,
                          unsigned int demand,
                          Verboseinfo *verboseinfo,
                          Error *err)
{
  unsigned int numofindexes;
  bool haserr = false;

  numofindexes = (unsigned int) strarray_size(indexnametab);
  emmesa->buf.nextaccessidx = emmesa->buf.nextstoreidx = 0;
  emmesa->numofindexes = numofindexes;
  emmesa->numofentries = numofindexes;
  emmesa->ident = (uint64_t) numofindexes;
  emmesa->trierep.encseqreadinfo = NULL;
  ALLOCASSIGNSPACE(emmesa->suffixarraytable,NULL,Suffixarray,numofindexes);
  ALLOCASSIGNSPACE(emmesa->nextpostable,NULL,Seqpos,numofindexes);
  if (inputthesequences(&emmesa->alpha,
                       emmesa->nextpostable,
                       emmesa->suffixarraytable,
                       indexnametab,
                       demand,
                       verboseinfo,
                       err) != 0)
  {
    haserr = true;
    return -1;
  }
  if (!haserr && numofindexes > 1U)
  {
    unsigned int idx;

    ALLOCASSIGNSPACE(emmesa->trierep.encseqreadinfo,NULL,Encseqreadinfo,
                     numofindexes);
    for (idx = 0; idx < numofindexes; idx++)
    {
      emmesa->trierep.encseqreadinfo[idx].encseqptr
        = emmesa->suffixarraytable[idx].encseq;
      emmesa->trierep.encseqreadinfo[idx].readmode
        = emmesa->suffixarraytable[idx].readmode;
    }
    initmergertrienodetable(&emmesa->trierep,(Seqpos) numofindexes,
                            numofindexes);
    if (insertfirstsuffixes(&emmesa->trierep,
                           emmesa->nextpostable,
                           emmesa->suffixarraytable,
                           numofindexes,
                           err) != 0)
    {
      FREESPACE(emmesa->trierep.encseqreadinfo);
      haserr = true;
    }
  }
  if (haserr)
  {
    FREESPACE(emmesa->suffixarraytable);
    FREESPACE(emmesa->nextpostable);
  }
  return haserr ? -1 : 0;
}

void wraptEmissionmergedesa(Emissionmergedesa *emmesa)
{
  unsigned int idx;

  for (idx = 0; idx < emmesa->numofindexes; idx++)
  {
    freesuffixarray(emmesa->suffixarraytable + idx);
  }
  FREESPACE(emmesa->suffixarraytable);
  FREESPACE(emmesa->trierep.encseqreadinfo);
  if (emmesa->numofindexes > 1U)
  {
    freemergertrierep(&emmesa->trierep);
  }
  FREESPACE(emmesa->nextpostable);
}
