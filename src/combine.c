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

#include "summarizer.h"

#if 0
static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static long get_long(const char * line, long * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%ld%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_string(const char * line, char ** value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");


  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except whitespace star and hash */
  char * end = start + strcspn(start, " \t\r\n*#");

  if (start == end) return 0;

  *end = 0;

  *value = xstrdup(start);
  free(s);

  return ws + end - start;
}

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long count_columns(const char * line)
{
  long columns = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tmp;
  size_t count = 0;
  size_t ws;


  while (1)
  {
    count = get_string(p, &tmp);
    if (!count) break;

    p += count;

    columns++;
    free(tmp);

    ws = strspn(p, " \t\r\n");
    if (!ws) break;

    p += ws;
  }

  free(s);
  return columns;
}
#endif

void cmd_combine()
{
  long i;
  long first = 1;
  long col_count = 0;
  long count;
  long sample_num;
  long line_count = 1;
  long file_line_count = 0;
  char * line;
  FILE * fp_list;
  FILE * fp_in;
  FILE * fp_out;
  FILE * fp_index;

  if (opt_skipcount != 1)
    fatal("Option --skip must be set to 1 when --combine");

  fp_list = xopen(opt_combine,"r");

  if (!opt_output)
    fatal("Option --combine requires an output file via --output");

  fp_out = xopen(opt_output,"w");

  char * s = NULL;
  xasprintf(&s, "%s.index", opt_output);
  fp_index = xopen(s,"w");
  free(s);

  while ((line=getnextline(fp_list)))
  {
    file_line_count = 0;
    char * filename = xstrdup(line);
    for (i = 0; i < strlen(filename); ++i)
      if (filename[i] == '\r' || filename[i] == '\n')
        filename[i] = 0;
    fp_in = xopen(filename,"r");

    fprintf(stdout, "Processing file %s\n", filename);

    line = getnextline(fp_in);

    /* if it was the first file */
    if (first)
    {
      col_count = count_columns(line);

      /* substract generations */
      col_count--;

      fprintf(fp_out,"%s\n",line);
    }
    else
    {
      long temp_col_count = count_columns(line);
      if (temp_col_count-1 != col_count)
        fatal("File %s contains %ld columns instead of %ld\n",
              temp_col_count, col_count+1);
    }

    /* we finished with the first file */
    if (first) first = 0;

    while ((line=getnextline(fp_in)))
    {
      double x;
      char * p = line;

      /* skip sample number */
      count = get_long(p,&sample_num);
      if (!count) goto l_unwind;

      p += count;

      fprintf(fp_out, "%ld", line_count);
      /* read remaining elements of current row */

      for (i = 0; i < col_count; ++i)
      {
        int decplaces = 0;
        count = get_double(p,&x, &decplaces);
        if (!count) goto l_unwind;

        p += count;

        fprintf(fp_out, "\t%.*f", decplaces,x);
      }

      line_count++;
      file_line_count++;
      fprintf(fp_out, "\n");
    }

    free(filename);
    fclose(fp_in);

    fprintf(fp_index,"%ld\n",file_line_count);
  }
l_unwind:
  fclose(fp_list);    
  fclose(fp_out);
  fclose(fp_index);
}
