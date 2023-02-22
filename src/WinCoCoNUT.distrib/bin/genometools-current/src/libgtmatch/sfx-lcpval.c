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
#include "seqpos-def.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "sfx-lcpval.h"

#include "sfx-cmpsuf.pr"
#define PRIu32
 struct Lcpvalueiterator
{
  Seqpos relpos,
         lastsuftabentry;
  Readmode readmode;
  const Encodedsequence *encseq;
  Encodedsequencescanstate *esr1, *esr2;
};

Lcpvalueiterator *newLcpvalueiterator(const Encodedsequence *encseq,
                                      Readmode readmode)
{
  Lcpvalueiterator *lvi;

  ALLOCASSIGNSPACE(lvi,NULL,Lcpvalueiterator,1);
  lvi->esr1 = newEncodedsequencescanstate();
  lvi->esr2 = newEncodedsequencescanstate();
  lvi->encseq = encseq;
  lvi->relpos = 0;
  lvi->readmode = readmode;
  lvi->lastsuftabentry = 0;
  return lvi;
}

Seqpos nextLcpvalueiterator(Lcpvalueiterator *lvi,
                            bool firstpage,
                            const Seqpos *suftabptr,
                            Seqpos numberofsuffixes)
{
  Seqpos lcpvalue;

  assert(lvi->relpos < numberofsuffixes);
  if (firstpage && lvi->relpos == 0)
  {
    lcpvalue = 0;
  } else
  {
    int cmp;

    cmp = comparetwosuffixes(lvi->encseq,
                             lvi->readmode,
                             &lcpvalue,
                             false,
                             false,
                             0,
                             lvi->lastsuftabentry,
                             suftabptr[lvi->relpos],
                             lvi->esr1,
                             lvi->esr2);
    if (cmp > 0)
    {
      fprintf(stderr,"pos=" FormatSeqpos
              ": cmp " FormatSeqpos
              " " FormatSeqpos " = %d, lcpval=" FormatSeqpos "\n",
              PRINTSeqposcast(lvi->relpos),
              PRINTSeqposcast(lvi->lastsuftabentry),
              PRINTSeqposcast(suftabptr[lvi->relpos]),
              cmp,
              PRINTSeqposcast(lcpvalue));
      exit(EXIT_FAILURE);
    }
  }
  lvi->lastsuftabentry = suftabptr[lvi->relpos];
  if (lvi->relpos + 1 == numberofsuffixes)
  {
    lvi->relpos = 0;
  } else
  {
    lvi->relpos++;
  }
  return lcpvalue;
}

void freeLcpvalueiterator(Lcpvalueiterator **lvi)
{
  freeEncodedsequencescanstate(&(*lvi)->esr1);
  freeEncodedsequencescanstate(&(*lvi)->esr2);
  FREESPACE(*lvi);
}
