/* -*- C++ -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library
 *
 * Copyright (C) 2003-2008
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_BITS_MINGW32_TIME_H
#define LEMON_BITS_MINGW32_TIME_H

#ifdef WIN32

#include <windows.h>
#include <ctime>
#include "dos.h"

char *asctime_r(const struct tm *t, char *buf);
struct tm * localtime_r (const time_t *t, struct tm *tm);
char *ctime_r(const time_t * tim_p , char * result);
int gettimeofday(struct timeval * tp, struct timezone *);

struct tms {
  long	tms_utime;
  long	tms_stime;
  long	tms_cutime;
  long	tms_cstime;
};

long filetime_to_clock(FILETIME *ft);

int times(struct tms *tmbuf);

#define _SC_CLK_TCK 1

int sysconf(int);

#endif

#endif
