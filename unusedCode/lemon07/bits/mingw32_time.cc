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

#ifdef WIN32

#include <lemon/bits/mingw32_time.h>

#include <windows.h>
#include <ctime>
#include "dos.h"

static const char days[] = 
"Sun Mon Tue Wed Thu Fri Sat ";
static const char months[] = 
"Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec ";

void num2str(char *c,int i) {
  c[0]=i/10+'0';
  c[1]=i%10+'0';
}

char *asctime_r(const struct tm *t, char *buf) {
  *(int*)buf=*(int*)(days+(t->tm_wday<<2));
  *(int*)(buf+4)=*(int*)(months+(t->tm_mon<<2));
  num2str(buf+8,t->tm_mday);
  if (buf[8]=='0') buf[8]=' ';
  buf[10]=' ';
  num2str(buf+11,t->tm_hour);
  buf[13]=':';
  num2str(buf+14,t->tm_min);
  buf[16]=':';
  num2str(buf+17,t->tm_sec);
  buf[19]=' ';
  num2str(buf+20,(t->tm_year+1900)/100);
  num2str(buf+22,(t->tm_year+1900)%100);
  buf[24]='\n'; buf[25]='\0';
  return buf;
}

struct tm * localtime_r (const time_t *t, struct tm *tm) {
  struct tm *tmp;
  
  if ((tmp = localtime(t)) && tm)
    *tm = *tmp;
  else
    return 0;
  
  return tm;
}

char *ctime_r(const time_t * tim_p , char * result) {
  struct tm tm;
  return asctime_r (localtime_r (tim_p, &tm), result);
}


int gettimeofday(struct timeval * tp, struct timezone *) {
  SYSTEMTIME systime;

  if (tp) {
    struct tm tmrec;
    time_t theTime = time(NULL);
    
    
    tmrec = *localtime(&theTime);
    tp->tv_sec = mktime(&tmrec);
    GetLocalTime(&systime); /* system time */

    tp->tv_usec = systime.wMilliseconds * 1000;
  }
  return 0;
}

long filetime_to_clock(FILETIME *ft)
{
  __int64 qw = ft->dwHighDateTime;
  qw <<= 32;
  qw |= ft->dwLowDateTime;
  qw /= 10000;  
  return (long) qw;

}

int times(struct tms *tmbuf)
{
  FILETIME create, exit, kernel, user;
  if (GetProcessTimes(GetCurrentProcess(),&create, &exit, &kernel, &user)) {
    tmbuf->tms_utime = filetime_to_clock(&user);
    tmbuf->tms_stime = filetime_to_clock(&kernel);
    tmbuf->tms_cutime = 0;
    tmbuf->tms_cstime = 0;
  }
  else {
    tmbuf->tms_utime = clock();
    tmbuf->tms_stime = 0;
    tmbuf->tms_cutime = 0;
    tmbuf->tms_cstime = 0;
  }
  return 0;
}

int sysconf(int)
{
  return 1;
}

#endif
