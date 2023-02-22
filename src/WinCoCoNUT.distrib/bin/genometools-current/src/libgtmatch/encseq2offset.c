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
#include "spacedef.h"
#include "sarr-def.h"
#include "safecast-gen.h"
#define PRIu32

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,Seqpos,Seqpos)

Seqpos *encseqtable2seqoffsets(Seqpos *totallength,
                               Specialcharinfo *specialcharinfo,
                               const Suffixarray *suffixarraytable,
                               unsigned int numofindexes)
{
  unsigned int idx;
  Uchar lastofprevious, firstofcurrent;
  Seqpos tmplength, *sequenceoffsettable;
  uint64_t tmpspecialcharacters,
           tmpspecialranges,
           tmprealspecialranges,
           tmplarge;

  assert(numofindexes > 0);
  ALLOCASSIGNSPACE(sequenceoffsettable,NULL,Seqpos,numofindexes);
  tmpspecialcharacters = (uint64_t) (numofindexes-1);
  tmpspecialranges = 0;
  tmprealspecialranges = 0;
  for (idx=0; idx<numofindexes; idx++)
  {
    if (idx == 0)
    {
      tmplength = 0;
      sequenceoffsettable[idx] = 0;
    } else
    {
      tmplength = getencseqtotallength(suffixarraytable[idx - 1].encseq);
      sequenceoffsettable[idx]
        = sequenceoffsettable[idx-1] + tmplength + (Seqpos) 1;
    }
    tmpspecialcharacters
      += (uint64_t) suffixarraytable[idx].specialcharinfo.specialcharacters;
    tmpspecialranges
      += (uint64_t) suffixarraytable[idx].specialcharinfo.specialranges;
    tmprealspecialranges
      += (uint64_t) suffixarraytable[idx].specialcharinfo.realspecialranges;
    if (idx > 0)
    {
      lastofprevious
        = getencodedchar(suffixarraytable[idx - 1].encseq, /* Random access */
                         tmplength-1,
                         suffixarraytable[idx - 1].readmode);
      firstofcurrent
        = getencodedchar(suffixarraytable[idx].encseq, /* Random access */
                         0,
                         suffixarraytable[idx].readmode);
      if (ISSPECIAL(lastofprevious))
      {
         if (ISSPECIAL(firstofcurrent))
         {
           tmpspecialranges--;
           tmprealspecialranges--;
         }
      } else
      {
        if (ISNOTSPECIAL(firstofcurrent))
        {
          tmpspecialranges++;
          tmprealspecialranges++;
        }
      }
    }
    tmplarge = (uint64_t) sequenceoffsettable[idx] +
               (uint64_t) getencseqtotallength(suffixarraytable[idx].encseq);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmplarge);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmpspecialcharacters);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmpspecialranges);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmprealspecialranges);
    printf("# seqlen[%u] = " FormatSeqpos "\n",
           idx,
           PRINTSeqposcast(getencseqtotallength(suffixarraytable[idx].encseq)));
  }
  tmplength = getencseqtotallength(suffixarraytable[numofindexes -1].encseq);
  *totallength = sequenceoffsettable[numofindexes-1] + tmplength;
  specialcharinfo->specialcharacters = (Seqpos) tmpspecialcharacters;
  specialcharinfo->specialranges = (Seqpos) tmpspecialranges;
  specialcharinfo->realspecialranges = (Seqpos) tmprealspecialranges;
  specialcharinfo->lengthofspecialprefix
    = suffixarraytable[0].specialcharinfo.lengthofspecialprefix;
  specialcharinfo->lengthofspecialsuffix
    = suffixarraytable[idx-1].specialcharinfo.lengthofspecialsuffix;
  return sequenceoffsettable;
}
