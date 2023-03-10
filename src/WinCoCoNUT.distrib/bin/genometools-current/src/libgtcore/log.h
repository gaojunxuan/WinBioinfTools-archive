/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef LOG_H
#define LOG_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

/* enable logging */
void log_enable(void);

/* returns true if logging is enabled, false otherwise */
bool log_enabled(void);

/* Prints the log message obtained from format and following parameters
   according if logging is enabled. The logging output is prefixed with the
   string "debug: " and finished by a newline.  */
void  log_log(const char *format, ...)
  __attribute__ ((format (printf, 1, 2)));

/* Prints the log message obtained from format and following parameter according
   to if logging is enabled analog to log_log(). But in contrast to
   log_log() log_vlog() does not accept individual arguments but a single
   va_list argument instead. */
void  log_vlog(const char *format, va_list);

/* return logging file pointer */
FILE* log_fp(void);

#endif
