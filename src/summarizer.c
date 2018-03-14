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

static char * progname;
static char progheader[80];
char * cmdline;

/* options */
long opt_help;
long opt_quiet;
long opt_skipcount;
long opt_version;
char * opt_combine;
char * opt_indexfile;
char * opt_output;
char * opt_summarize;

static struct option long_options[] =
{
  {"help",       no_argument,       0, 0 },  /* 0 */
  {"version",    no_argument,       0, 0 },  /* 1 */
  {"quiet",      no_argument,       0, 0 },  /* 2 */
  {"summarize",  required_argument, 0, 0 },  /* 3 */
  {"combine",    required_argument, 0, 0 },  /* 4 */
  {"output",     required_argument, 0, 0 },  /* 5 */
  {"skip",       required_argument, 0, 0 },  /* 6 */
  {"index",      required_argument, 0, 0 },  /* 7 */
  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;

  /* set defaults */

  progname = argv[0];

  opt_summarize = NULL;
  opt_combine = NULL;
  opt_indexfile = NULL;
  opt_output = NULL;
  opt_help = 0;
  opt_version = 0;
  opt_quiet = 0;
  opt_skipcount = 1;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        opt_summarize = xstrdup(optarg);
        break;

      case 4:
        opt_combine = xstrdup(optarg);
        break;

      case 5:
        opt_output = xstrdup(optarg);
        break;

      case 6:
        opt_skipcount = atol(optarg);
        if (opt_skipcount < 0)
          fatal("option --skip requires a positive integer or 0");
        break;

      case 7:
        opt_indexfile = xstrdup(optarg);
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands = 0;

  /* check for number of independent commands selected */

  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_summarize)
    commands++;
  if (opt_combine)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  if (!commands)
  {
    opt_help = 1;
    return;
  }
}

static void dealloc_switches()
{
  if (opt_summarize) free(opt_summarize);
  if (opt_combine) free(opt_combine);
  if (opt_indexfile) free(opt_indexfile);
  if (opt_output) free(opt_output);
}

void cmd_help()
{
  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help                display help information\n"
          "  --version             display version information\n"
          "  --quiet               only output warnings and fatal errors to stderr\n"
          "  --summarize FILENAME  summarize MCMC file\n"
          "  --combine FILENAME    combine list of MCMC files in specified file\n"
          "  --output FILENAME     write output to specified file\n"
          "  --skip INTEGER        skip INTEGER lines from beginning of MCMC files\n"
          "\n"
         );

  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           arch_get_cores());
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/summarizer\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_version)
  {
    ;
  }
  else if (opt_summarize)
  {
    if (opt_indexfile)
      cmd_summary_full();
    else
      cmd_summary();
  }
  else if (opt_combine)
  {
    cmd_combine();
  }

  dealloc_switches();
  free(cmdline);
  return (0);
}
