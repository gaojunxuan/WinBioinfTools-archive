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
#include <errno.h>
#include <string.h>
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/str.h"
#include "ushort-def.h"
#include "fmi-bwtbound.h"
#include "intbits.h"
#include "safecast-gen.h"
#include "mapspec-def.h"
#include "format64.h"
#include "stamp.h"
#define PRIu32
#define PRIu64
#define ASSIGNPTR2STARTPTR(TYPE)\
        if (mapspec->numofunits == 0)\
        {\
          *((TYPE **) mapspec->startptr) = NULL;\
        } else\
        {\
          voidptr = (((char *) ptr) + byteoffset);\
          *((TYPE **) mapspec->startptr) = voidptr;\
        }

#define ALIGNSIZE sizeof (void *)

#define WRITEACTIONWITHTYPE(TYPE)\
        if (fwrite(*((TYPE **) mapspecptr->startptr),\
                   mapspecptr->sizeofunit,\
                   (size_t) mapspecptr->numofunits, fp) !=\
                    (size_t) mapspecptr->numofunits)\
        {\
          error_set(err,"cannot write %lu items of size %u: "\
                            "errormsg=\"%s\"",\
                        (unsigned long) mapspecptr->numofunits,\
                        (unsigned int) mapspecptr->sizeofunit,\
                        strerror(errno));\
          haserr = true;\
        }

static uint64_t detexpectedaccordingtomapspec(const ArrayMapspecification
                                              *mapspectable)
{
  uint64_t sumup = 0;
  Mapspecification *mapspecptr;

  for (mapspecptr = mapspectable->spaceMapspecification;
       mapspecptr < mapspectable->spaceMapspecification +
                    mapspectable->nextfreeMapspecification; mapspecptr++)
  {
    sumup += (uint64_t) mapspecptr->sizeofunit *
             (uint64_t) mapspecptr->numofunits;
    if (sumup % ALIGNSIZE > 0)
    {
      sumup += (ALIGNSIZE - (sumup % ALIGNSIZE));
    }
  }
  return sumup;
}

#undef SKDEBUG
#ifdef SKDEBUG
static void showmapspec(const Mapspecification *mapspec)
{
  printf("(%s,size=%lu,elems=%lu)",
           mapspec->name,
           (unsigned long) mapspec->sizeofunit,
           mapspec->numofunits);
}
#endif

static int assigncorrecttype(Mapspecification *mapspec,
                             void *ptr,
                             unsigned long byteoffset,
                             Error *err)
{
  void *voidptr;
  bool haserr = false;

  error_check(err);
  switch (mapspec->typespec)
  {
    case UcharType:
      ASSIGNPTR2STARTPTR(Uchar);
      break;
    case UshortType:
      ASSIGNPTR2STARTPTR(Ushort);
      break;
    case Uint32Type:
      ASSIGNPTR2STARTPTR(uint32_t);
      break;
    case UnsignedlongType:
      ASSIGNPTR2STARTPTR(Unsignedlong);
      break;
    case Uint64Type:
      ASSIGNPTR2STARTPTR(uint64_t);
      break;
    case BitstringType:
      ASSIGNPTR2STARTPTR(Bitstring);
      break;
    case SeqposType:
      ASSIGNPTR2STARTPTR(Seqpos);
      break;
    case BwtboundType:
      ASSIGNPTR2STARTPTR(Bwtbound);
      break;
    case PairBwtidxType:
      ASSIGNPTR2STARTPTR(PairBwtidx);
      break;
    case TwobitencodingType:
      ASSIGNPTR2STARTPTR(Twobitencoding);
      break;
    default:
      error_set(err,"no assignment specification for size %lu",
                    (unsigned long) mapspec->sizeofunit);
      haserr = true;
  }
  return haserr ? -1 : 0;
}

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,unsigned long,unsigned_long)

int fillmapspecstartptr(Assignmapspec assignmapspec,
                        void **mappeduserptr,
                        void *assignmapinfo,
                        const Str *tmpfilename,
                        unsigned long expectedsize,
                        Error *err)
{
  void *mapptr;
  uint64_t expectedaccordingtomapspec;
  unsigned long byteoffset = 0;
  size_t numofbytes;
  ArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  bool haserr = false;
  unsigned long totalpadunits = 0;

  error_check(err);
  INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,false);
  mapptr = fa_mmap_read(str_get(tmpfilename), &numofbytes);
  if (mapptr == NULL)
  {
    error_set(err,"could not map datafile %s",str_get(tmpfilename));
    haserr = true;
  }
  *mappeduserptr = mapptr;
  if (!haserr)
  {
    if (assigncorrecttype(mapspectable.spaceMapspecification,
                          mapptr,0,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    expectedaccordingtomapspec = detexpectedaccordingtomapspec(&mapspectable);
    if (expectedaccordingtomapspec != (uint64_t) numofbytes)
    {
      error_set(err,"%lu bytes read from %s, but " Formatuint64_t
                         " expected",
                         (unsigned long) numofbytes,
                         str_get(tmpfilename),
                         PRINTuint64_tcast(expectedaccordingtomapspec));
      haserr = true;
    }
  }
  if (!haserr)
  {
    mapspecptr = mapspectable.spaceMapspecification;
    assert(mapspecptr != NULL);
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (byteoffset % (unsigned long) ALIGNSIZE > 0)
    {
      size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
      byteoffset += (unsigned long) padunits;
      totalpadunits += (unsigned long) padunits;
    }
    for (mapspecptr++;
         mapspecptr < mapspectable.spaceMapspecification +
                      mapspectable.nextfreeMapspecification; mapspecptr++)
    {
      if (assigncorrecttype(mapspecptr,mapptr,byteoffset,err) != 0)
      {
        haserr = true;
        break;
      }
      byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                                (uint64_t) (byteoffset +
                                            mapspecptr->sizeofunit *
                                            mapspecptr->numofunits));
      if (byteoffset % (unsigned long) ALIGNSIZE > 0)
      {
        size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
        byteoffset += (unsigned long) padunits;
        totalpadunits += (unsigned long) padunits;
      }
    }
  }
  if (!haserr)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      error_set(err,"mapping: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        expectedsize,byteoffset);
      haserr = true;
    }
  }
  FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}

int flushtheindex2file(FILE *fp,
                       Assignmapspec assignmapspec,
                       void *assignmapinfo,
                       unsigned long expectedsize,
                       Error *err)
{
  ArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  unsigned long byteoffset = 0;
  bool haserr = false;
  Uchar padbuffer[ALIGNSIZE-1] = {0};
  unsigned long totalpadunits = 0;

  error_check(err);
  INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,true);
  assert(mapspectable.spaceMapspecification != NULL);
  for (mapspecptr = mapspectable.spaceMapspecification;
       mapspecptr < mapspectable.spaceMapspecification +
                    mapspectable.nextfreeMapspecification;
       mapspecptr++)
  {
#ifdef SKDEBUG
    printf("# flushtheindex2file");
    showmapspec(mapspecptr);
    printf(" at byteoffset %lu\n",byteoffset);
#endif
    if (mapspecptr->numofunits > 0)
    {
      switch (mapspecptr->typespec)
      {
        case UcharType:
          WRITEACTIONWITHTYPE(Uchar);
          break;
        case UshortType:
          WRITEACTIONWITHTYPE(Ushort);
          break;
        case Uint32Type:
          WRITEACTIONWITHTYPE(uint32_t);
          break;
        case UnsignedlongType:
          WRITEACTIONWITHTYPE(Unsignedlong);
          break;
        case Uint64Type:
          WRITEACTIONWITHTYPE(uint64_t);
          break;
        case BitstringType:
          WRITEACTIONWITHTYPE(Bitstring);
          break;
        case SeqposType:
          WRITEACTIONWITHTYPE(Seqpos);
          break;
        case BwtboundType:
          WRITEACTIONWITHTYPE(Bwtbound);
          break;
        case PairBwtidxType:
          WRITEACTIONWITHTYPE(PairBwtidx);
          break;
        case TwobitencodingType:
          WRITEACTIONWITHTYPE(Twobitencoding);
          break;
        default:
           error_set(err,"no map specification for size %lu",
                         (unsigned long) mapspecptr->sizeofunit);
           haserr = true;
      }
    }
    if (haserr)
    {
      break;
    }
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (byteoffset +
                                          mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (byteoffset % (unsigned long) ALIGNSIZE > 0)
    {
      size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
      if (fwrite(padbuffer,
                sizeof (Uchar),padunits,fp) != padunits)
      {
        error_set(err,"cannot write %lu items of size %u: "
                          "errormsg=\"%s\"",
                           (unsigned long) padunits,
                           (unsigned int) sizeof (Uchar),
                           strerror(errno));
        haserr = true;
      }
      byteoffset += (unsigned long) padunits;
      totalpadunits += (unsigned long) padunits;
    }
  }
  if (!haserr)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      error_set(err,"flushindex: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        expectedsize,
                        byteoffset);
      haserr = true;
    }
  }
  FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}
