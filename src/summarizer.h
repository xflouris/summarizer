/*
    Copyright (C) 2018 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, United Kingdom
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

typedef struct snode_s
{
  char * label;
  double length;
  struct snode_s * left;
  struct snode_s * right;
  struct snode_s * parent;
  unsigned int leaves;

  void * data;

  int mark;
  double age;

  unsigned int node_index;
  unsigned int diploid;

} snode_t;

typedef struct stree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  snode_t ** nodes;

  snode_t * root;
} stree_t;

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
extern long opt_map_median;
extern long opt_map_hpdci;
extern char * opt_combine;
extern char * opt_indexfile;
extern char * opt_mapfile;
extern char * opt_output;
extern char * opt_summarize;
extern char * opt_treefile;

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
long get_double(const char * line, double * value, int * decplaces);
long get_string(const char * line, char ** value);
void reallocline(size_t newmaxsize);
char * getnextline(FILE * fp);
long count_columns(const char * line);

/* functions in parse_stree.y */

void stree_destroy(stree_t * tree,
                   void (*cb_destroy)(void *));
stree_t * stree_parse_newick(const char * filename);
stree_t * stree_parse_newick_string(const char * s);

/* functions in stree.c */

char * stree_export_newick(const snode_t * root,
                           char * (*cb_serialize)(const snode_t *));

/* functions in map.c */

void cmd_map(void);
