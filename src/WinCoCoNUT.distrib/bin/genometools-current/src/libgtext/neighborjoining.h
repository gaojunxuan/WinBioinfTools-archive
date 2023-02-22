/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef NEIGHBORJOINING_H
#define NEIGHBORJOINING_H

#include <stdio.h>
#include "libgtcore/error.h"

typedef struct NeighborJoining NeighborJoining;

typedef double (*NeighborJoiningDistFunc)(unsigned long, unsigned long, void*);

NeighborJoining* neighborjoining_new(unsigned long num_of_taxa, void *data,
                                     NeighborJoiningDistFunc);
void             neighborjoining_show_tree(const NeighborJoining*, FILE*);
void             neighborjoining_delete(NeighborJoining*);

#endif
