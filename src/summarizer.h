/*
    copyright (c) 2018 tomas flouri

    this program is free software: you can redistribute it and/or modify
    it under the terms of the gnu affero general public license as
    published by the free software foundation, either version 3 of the
    license, or (at your option) any later version.

    this program is distributed in the hope that it will be useful,
    but without any warranty; without even the implied warranty of
    merchantability or fitness for a particular purpose.  see the
    gnu affero general public license for more details.

    you should have received a copy of the gnu affero general public license
    along with this program.  if not, see <http://www.gnu.org/licenses/>.

    contact: tomas flouri <t.flouris@ucl.ac.uk>,
    department of genetics, evolution and environment, 
    university college london, gower street, wc1e 6bt london, united kingdom
*/

#include <assert.h>
#include <fcntl.h>
#include <getopt.h>
#include <search.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/stat.h>
#include <stdint.h>

#ifndef _MSC_VER
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

/* constants */

#define PROG_NAME "summarizer"
#define PROG_VERSION "v0.0.1"

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#else
#error "Big endian ppc64 CPUs not supported"
#endif

#else

#define PROG_CPU "x86_64"

#endif

#ifdef __APPLE__
#define PROG_OS "osx"
#endif

#ifdef __linux__
#define PROG_OS "linux"
#endif

#ifdef _WIN32
#define PROG_OS "win"
#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU

#define LINEALLOC 2048

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef _MSC_VER
#define SWAP(x,y) do                                                  \
  {                                                                   \
    size_t s = MAX(sizeof(x),sizeof(y));                              \
    unsigned char * temp = (unsigned char *)malloc(s*sizeof(char));   \
    memcpy(temp,&y,s);                                                \
    memcpy(&y,&x,s);                                                  \
    memcpy(&x,temp,s);                                                \
    free(temp);                                                       \
  } while(0)
#else
#define SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)
#endif

/* options */

extern long opt_help;
extern long opt_quiet;
extern long opt_skipcount;
extern long opt_version;
extern char * opt_combine;
extern char * opt_indexfile;
extern char * opt_output;
extern char * opt_summarize;

/* functions in summary.c */

void cmd_summary(void);

/* functions in summaryfull.c */

void cmd_summary_full(void);

/* functions in combine.c */

void cmd_combine(void);

/* functions in util.c */

#ifdef _MSC_VER
int xasprintf(char ** strp, const char * fmt, ...);
__declspec(noreturn) void fatal(const char * format, ...);
#else
void fatal(const char * format, ...) __attribute__ ((noreturn));
#define xasprintf asprintf
#endif
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
FILE * xopen(const char * filename, const char * mode);
void * pll_aligned_alloc(size_t size, size_t alignment);
void pll_aligned_free(void * ptr);

/* functions in arch.c */

unsigned long arch_get_memused(void);

unsigned long arch_get_memtotal(void);

long arch_get_cores(void);

/* functions in parse.c */

long get_long(const char * line, long * value);
long get_double(const char * line, double * value);
long get_string(const char * line, char ** value);
void reallocline(size_t newmaxsize);
char * getnextline(FILE * fp);
long count_columns(const char * line);
