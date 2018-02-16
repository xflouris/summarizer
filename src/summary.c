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

void cmd_summary()
{
  FILE * fp;
  long i,j;
  long count;
  long opt_samples;
  long sample_num;
  long line_count = 0;
  char * line;
  FILE * fp_out;

  if (opt_skipcount < 1)
    fatal("Option --skip must be greater or equal to 1");

  opt_samples = count_lines(opt_summarize);

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
      count = get_double(p,&x);
      if (!count) goto l_unwind;

      p += count;

      matrix[i][line_count] = x;
    }

    line_count++;
  }
  progress_done();

  fprintf(fp_out, "%s",labels[1]);
  for (i = 1; i < col_count; ++i)
    fprintf(fp_out, " %s", labels[i+1]);
  fprintf(fp_out,"\n");

  /* compute means */
  printf("Computing means...\n");
  fprintf(fp_out, "mean    ");
  for (i = 0; i < col_count; ++i)
  {
    double sum = 0;
    for (j = 0; j < opt_samples; ++j)
      sum += matrix[i][j];

    mean[i] = sum/(opt_samples);
    fprintf(fp_out, "  %f", mean[i]);
  }
  fprintf(fp_out, "\n");

  /* compute standard deviation */
  printf("Computing stdev...\n");
  for (i = 0; i < col_count; ++i)
  {
    double sd = 0;
    for (j = 0; j < opt_samples; ++j)
      sd += (matrix[i][j]-mean[i]) * (matrix[i][j]-mean[i]);

    stdev[i] = sqrt(sd/(opt_samples-1));
  }

  #ifdef COMPUTE_ESS
  /* compute tint */
  printf("Computing tint...\n");
  for (i = 0; i < col_count; ++i)
    tint[i] = eff_ict(matrix[i],opt_samples,mean[i],stdev[i]);
  #endif

  /* compute and print medians */
  printf("Computing medians...\n");
  fprintf(fp_out, "median  ");
  long median_line = opt_samples / 2;

  for (i = 0; i < col_count; ++i)
  {
    qsort(matrix[i], opt_samples, sizeof(double), cb_cmp_double);

    double median = matrix[i][median_line];
    if ((opt_samples & 1) == 0)
    {
      median += matrix[i][median_line-1];
      median /= 2;
    }

    fprintf(fp_out, "  %f", median);
  }
  fprintf(fp_out, "\n");

  /* print standard deviation */
  if (opt_output)
    printf("Writing output to %s...\n", opt_output);
  else
    printf("Writing output...\n");

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
    fprintf(fp_out, "  %f", matrix[i][0]);
  }
  fprintf(fp_out, "\n");

  /* print maximum values */
  fprintf(fp_out, "max     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][opt_samples-1]);
  }
  fprintf(fp_out, "\n");

  /* print line at 2.5% of matrix */
  fprintf(fp_out, "2.5%%    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][(long)(opt_samples*.025)]);
  }
  fprintf(fp_out, "\n");

  /* print line at 97.5% of matrix */
  fprintf(fp_out, "97.5%%   ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "  %f", matrix[i][(long)(opt_samples*.975)]);
  }
  fprintf(fp_out, "\n");

  /* compute and print HPD 2.5% and 97.5% */
  for (i = 0; i < col_count; ++i)
    hpd_interval(matrix[i],opt_samples,hpd025+i,hpd975+i,0.05);

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
    fprintf(fp_out, "  %f", opt_samples/tint[i]);
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

  /* table-like summary */
  fprintf(fp_out, "\n\nPosterior mean (95%% Equal-tail CI) (95%% HPD CI) HPD-CI-width\n\n");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(fp_out, "%-15s ", labels[i+1]);
    fprintf(fp_out, "%f ",mean[i]);
    fprintf(fp_out, "(%f, %f) ", matrix[i][(long)(opt_samples*.025)],
                                 matrix[i][(long)(opt_samples*.975)]);
    fprintf(fp_out, "(%f, %f) ", hpd025[i], hpd975[i]);
    fprintf(fp_out, "%f\n", hpd975[i] - hpd025[i]);
  }
  for (i = 0; i < col_count+1; ++i)
    free(labels[i]);
  free(labels);

l_unwind:
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
