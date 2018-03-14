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

static char ** getlabels(const char * line)
{
  long columns = 0;
  size_t count = 0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;
  char * tmp;
  char ** labels;


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

  labels = (char **)xmalloc((size_t)columns * sizeof(char *));
  columns = 0;
  p = s;
  while (1)
  {
    count = get_string(p, &tmp);
    if (!count) break;

    p += count;

    labels[columns++] = tmp;

    ws = strspn(p, " \t\r\n");
    if (!ws) break;

    p += ws;
  }

  free(s);

  return labels;
}

static long countindexrecords(long * dataset_records_count)
{
  long count;
  long linelong = 0;
  long line_count = 0;
  long records = 0;
  char * line;
  char * p;
  FILE * fp;

  fp = xopen(opt_indexfile,"r");

  while((line=getnextline(fp)))
  {
    ++line_count;
    p = line;

    count = get_long(p,&linelong);
    if (!count)
      fatal("Invalid entry in line %ld of %s", line_count, opt_indexfile);

    if (linelong <= 0)
      fatal("Invalid record on line %ld of %s", line_count, opt_indexfile);
    dataset_records_count[line_count-1] = linelong;
    records += linelong;
  }

  fclose(fp);

  return records;
}

static long count_lines(const char * mcmcfile)
{
  FILE * fp;
  long line_count = 0;

  fp = xopen(mcmcfile,"r");

  while(getnextline(fp))
    ++line_count;

  fclose(fp);

  return line_count;
}

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  return *x > *y;
}

#ifdef COMPUTE_ESS
static double eff_ict(double * y, long n, double mean, double stdev)
{
  /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
     sequence method */

  long i,j;
  double tint = 1;
  double rho, rho0 = 0;

  /* TODO: ADDED NOW */
  double * x = (double *)xmalloc((size_t)n * sizeof(double));
  for (i = 0; i < n; ++i)
    x[i] = (y[i]-mean)/stdev;


  if (stdev/(fabs(mean)+1) < 1E-9)
  {
   tint = n;
  }
  else
  {
    for (i = 1; i < n-10; ++i)
    {
      rho = 0;
      for (j = 0; j < n - i; ++j)
        rho += x[j]*x[i+j];

      rho /= (n-1);

      if (i > 10 && rho+rho0 < 0)
        break;

      tint += rho*2;
      rho0 = rho;
    }
  }

  free(x);

  return tint;
}
#endif

static void hpd_interval(double * x,
                         long n,
                         double * ltail,
                         double * rtail,
                         double alpha)
{
  long lrow = (long)(n*alpha/2);
  long urow = (long)(n*(1-alpha/2));
  long diffrow = urow - lrow;
  long l,r;

  long left = lrow;
  

  double w = x[urow] - x[lrow];

  *ltail = x[lrow];
  *rtail = x[urow];

  if (n <= 2) return;

  for (l=0,r=l+diffrow; r < n; l++,r++)
  {
    if (x[r] - x[l] < w)
    {
      left = l;
      w = x[r] - x[l];
    }
  }

  *ltail = x[left];
  *rtail = x[left + diffrow];
}

static void print_summary(long index,
                          long start,
                          long records,
                          long col_count,
                          char ** labels,
                          double ** matrix)
{
  long i,j;
  FILE * fp_out;

  double * mean = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * medianarray = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd025 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)col_count * sizeof(double));
  #ifdef COMPUTE_ESS
  double * tint = (double *)xmalloc((size_t)col_count * sizeof(double));
  #endif
  double * stdev = (double *)xmalloc((size_t)col_count * sizeof(double));

  char * s = NULL;

  if (index)
    xasprintf(&s, "%s.%ld.txt", opt_output, index);
  else
    xasprintf(&s, "%s.combined.txt", opt_output);

  fp_out = xopen(s,"w");
  if (index)
    printf("Summarizing dataset %ld in %s\n", index, s);
  else
    printf("Summarizing combined dataset in %s\n", s);
  free(s);


  fprintf(fp_out, "%s",labels[1]);
  for (i = 1; i < col_count; ++i)
    fprintf(fp_out, " %s", labels[i+1]);
  fprintf(fp_out,"\n");

  /* compute means */
  #if 0
  printf("Computing means...\n");
  #endif
  fprintf(fp_out, "mean    ");
  for (i = 0; i < col_count; ++i)
  {
    double sum = 0;
    for (j = start; j < start+records; ++j)
      sum += matrix[i][j];

    mean[i] = sum/records;
    fprintf(fp_out, "  %f", mean[i]);
  }
  fprintf(fp_out, "\n");

  /* compute standard deviation */
  #if 0
  printf("Computing stdev...\n");
  #endif
  for (i = 0; i < col_count; ++i)
  {
    double sd = 0;
    for (j = start; j < start+records; ++j)
      sd += (matrix[i][j]-mean[i]) * (matrix[i][j]-mean[i]);

    stdev[i] = sqrt(sd/(records-1));
  }

  #ifdef COMPUTE_ESS
  /* compute tint */
  #if 0
  printf("Computing tint...\n");
  #endif
  for (i = 0; i < col_count; ++i)
    tint[i] = eff_ict(matrix[i]+start,records,mean[i],stdev[i]);
  #endif

  /* compute and print medians */
  #if 0
  printf("Computing medians...\n");
  #endif
  fprintf(fp_out, "median  ");
  long median_line = records / 2;

  for (i = 0; i < col_count; ++i)
  {
    qsort(matrix[i]+start, records, sizeof(double), cb_cmp_double);

    double median = matrix[i][start+median_line];
    if ((records & 1) == 0)
    {
      median += matrix[i][start+median_line-1];
      median /= 2;
    }

    fprintf(fp_out, "  %f", median);
    medianarray[i] = median;
  }
  fprintf(fp_out, "\n");

  /* print standard deviation */
  #if 0
  if (opt_output)
    printf("Writing output to %s...\n", opt_output);
  else
    printf("Writing output...\n");
  #endif

  fprintf(fp_out, "S.D     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", stdev[i]);
  }
  fprintf(fp_out, "\n");

  /* print minimum values */
  fprintf(fp_out, "min     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][start]);
  }
  fprintf(fp_out, "\n");

  /* print maximum values */
  fprintf(fp_out, "max     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][start+records-1]);
  }
  fprintf(fp_out, "\n");

  /* print line at 2.5% of matrix */
  fprintf(fp_out, "2.5%%    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][start+(long)(records*.025)]);
  }
  fprintf(fp_out, "\n");

  /* print line at 97.5% of matrix */
  fprintf(fp_out, "97.5%%   ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][start+(long)(records*.975)]);
  }
  fprintf(fp_out, "\n");

  /* compute and print HPD 2.5% and 97.5% */
  for (i = 0; i < col_count; ++i)
    hpd_interval(matrix[i]+start,records,hpd025+i,hpd975+i,0.05);

  /* print 2.5% HPD */
  fprintf(fp_out, "2.5%%HPD ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", hpd025[i]);
  }
  fprintf(fp_out, "\n");

  /* print 97.5% HPD */
  fprintf(fp_out, "97.5%%HPD");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", hpd975[i]);
  }
  fprintf(fp_out, "\n");

  #ifdef COMPUTE_ESS
  /* print ESS */
  fprintf(fp_out, "ESS*    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", records/tint[i]);
  }
  fprintf(fp_out, "\n");
    
  /* print Eff */
  fprintf(fp_out, "Eff*    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", 1/tint[i]);
  }
  fprintf(fp_out, "\n");
  #endif

  fclose(fp_out);

  s = NULL;
  if (index)
    xasprintf(&s, "%s.table.%ld.txt", opt_output, index);
  else
    xasprintf(&s, "%s.table.combined.txt", opt_output);

  fp_out = xopen(s,"w");
  printf("Table-like summary in %s\n", s);
  free(s);

  /* table-like summary */
  fprintf(fp_out, "Posterior median mean (95%% Equal-tail CI) (95%% HPD CI) HPD-CI-width\n");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "%-15s ", labels[i+1]);
    fprintf(fp_out, "%f ",medianarray[i]);
    fprintf(fp_out, "%f ",mean[i]);
    fprintf(fp_out, "(%f, %f) ", matrix[i][start+(long)(records*.025)],
                                 matrix[i][start+(long)(records*.975)]);
    fprintf(fp_out, "(%f, %f) ", hpd025[i], hpd975[i]);
    fprintf(fp_out, "%f\n", hpd975[i] - hpd025[i]);
  }

  free(mean);
  free(medianarray);
  free(hpd025);
  free(hpd975);
  #ifdef COMPUTE_ESS
  free(tint);
  #endif
  free(stdev);

  fclose(fp_out);
}

void cmd_summary_full()
{
  FILE * fp;
  long i;
  long count;
  long opt_samples;
  long sample_num;
  long line_count = 0;
  long dataset_count = 0;
  char * line;
  FILE * fp_out;
  long * dataset_records_count = NULL;

  if (opt_skipcount < 1)
    fatal("Option --skip must be greater or equal to 1");

  opt_samples = count_lines(opt_summarize);
  dataset_count = count_lines(opt_indexfile);

  printf("Reading records...\n");

  dataset_records_count = (long *)xcalloc(dataset_count,sizeof(long));
  long total_records_count = countindexrecords(dataset_records_count);

  if (opt_samples-opt_skipcount != total_records_count)
    fatal("Number of records in %s does not match with index file %s",
          opt_summarize, opt_indexfile);

  /* skip header */
  opt_samples -= 1;

  fp = xopen(opt_summarize,"r");

  if (opt_output)
    fp_out = xopen(opt_output,"w");
  else
    fp_out = stdout;

  /* skip line containing header */
  line = getnextline(fp);

  /* compute number of columns in the file */
  long col_count = 0;

  col_count = count_columns(line);
  char ** labels = getlabels(line);

  /* subtract generations */
  col_count--;

  for (i = 1; i < opt_skipcount; ++i)
    line = getnextline(fp);


  fprintf(stdout, "Skipped %ld header line(s)...\n", opt_skipcount);
  fprintf(stdout, "Processing %ld lines (samples) each %ld columns...\n",
          opt_samples, col_count);

  double ** matrix = (double **)xmalloc((size_t)col_count * sizeof(double *));
  for (i = 0; i < col_count; ++i)
    matrix[i] = (double *)xmalloc((size_t)opt_samples * sizeof(double));

  double * mean = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd025 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)col_count * sizeof(double));
  #ifdef COMPUTE_ESS
  double * tint = (double *)xmalloc((size_t)col_count * sizeof(double));
  #endif
  double * stdev = (double *)xmalloc((size_t)col_count * sizeof(double));

  unsigned long total_steps = opt_samples;
  progress_init("Processing data...", total_steps);

  while((line=getnextline(fp)))
  {
    double x;
    char * p = line;

    progress_update(line_count);

    /* skip sample number */
    count = get_long(p,&sample_num);
    if (!count) goto l_unwind;

    p += count;

    /* read remaining elements of current row */

    for (i = 0; i < col_count; ++i)
    {
      count = get_double(p,&x,NULL);
      if (!count) goto l_unwind;

      p += count;

      matrix[i][line_count] = x;
    }

    line_count++;
  }
  progress_done();

  /* print individual dataset summaries */
  long start = 0;
  for (i = 0; i < dataset_count; ++i)
  {
    print_summary(i+1,
                  start,
                  dataset_records_count[i],
                  col_count,
                  labels,
                  matrix);

    start += dataset_records_count[i];
  }

  /* print combined summary */
  printf("Summarizing combined dataset...\n");
  print_summary(0,
                0,
                opt_samples,
                col_count,
                labels,
                matrix);

l_unwind:
  for (i = 0; i < col_count+1; ++i)
    free(labels[i]);
  free(labels);

  if (dataset_records_count)
    free(dataset_records_count);

  for (i = 0; i < col_count; ++i)
    free(matrix[i]);
  free(matrix);

  free(mean);
  free(hpd025);
  free(hpd975);
  #ifdef COMPUTE_ESS
  free(tint);
  #endif
  free(stdev);

  if (opt_output)
    fclose(fp_out);

  fclose(fp);
}
