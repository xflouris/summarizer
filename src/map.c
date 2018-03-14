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

static void setbranchlengths(stree_t * tree)
{
  long i;
  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
    if (!tree->nodes[i]->mark)
      fatal("No entries in %s found for some inner nodes (%s)",
            opt_mapfile,
            tree->nodes[i]->label ? tree->nodes[i]->label :
              "Some internal nodes in the tree file have no label");

  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    snode_t * node = tree->nodes[i];

    if (node->parent)
      node->length = node->parent->age - node->age;
    else
      node->length = 0;
  }
}

void cmd_map()
{
  long i;
  long count;
  char * line;
  char * label;
  FILE * fp_out;

  if (!opt_treefile)
    fatal("Missing tree file (option --tree)");

  stree_t * t = stree_parse_newick(opt_treefile);

  FILE * fp_table = xopen(opt_mapfile, "r");


  /* skip header line from table */
  line = getnextline(fp_table);

  for (i = 0; i < t->tip_count + t->inner_count; ++i)
  {
    t->nodes[i]->mark = 0;
    t->nodes[i]->age = 0;
  }

  if (opt_output)
    fp_out = xopen(opt_output, "w");
  else
    fp_out = stdout;

  while ((line=getnextline(fp_table)))
  {
    char * p = line;

    /* remove all commas, opening and closing parentheses */
    while (*p)
    {
      if (*p == '(' || *p == ')' || *p == ',')
        *p = ' ';
      ++p;
    }

    p = line;

    /* skip all white-space */
    size_t ws = strspn(p, " \t\r\n");

    /* is it a blank line or comment ? */
    if (!p[ws] || p[ws] == '#')
      continue;

    p += ws;

    count = get_string(p,&label);
    if (!count)
      fatal("File %s contains no data...", opt_mapfile);

    p += count;

    if ((strlen(label) <= 3) || strncmp(label,"t_n",3))
    {
      free(label);
      break;
    }

    for (i = t->tip_count; i < t->tip_count + t->inner_count; ++i)
    {
      if (!strcmp(t->nodes[i]->label,label+3))
      {
        break;
      }
    }
    if (i == t->tip_count + t->inner_count)
      fatal("Could not find label %s\n", label+3);

    snode_t * node = t->nodes[i];

    if (node->mark)
      fatal("Duplicate record for node %s", t->nodes[i]->label);

    double median;
    double mean;
    double etlo;
    double ethi;
    double hpdlo;
    double hpdhi;

    count = get_double(p,&median,NULL);
    if (!count)
      fatal("Cannot read median for node %s", label+3);

    p += count;

    count = get_double(p,&mean,NULL);
    if (!count)
      fatal("Cannot read mean for node %s", label+3);
    p += count;
    
    count = get_double(p,&etlo,NULL);
    if (!count)
      fatal("2.Cannot read equal-tail CI for node %s", label+3);
    p += count;

    count = get_double(p,&ethi,NULL);
    if (!count)
      fatal("4.Cannot read equal-tail CI for node %s", label+3);
    p += count;

    count = get_double(p,&hpdlo,NULL);
    if (!count)
      fatal("Cannot read HPD CI for node %s", label+3);
    p += count;

    count = get_double(p,&hpdhi,NULL);
    if (!count)
      fatal("Cannot read HPD CI for node %s", label+3);
    p += count;


    #if 0
    printf("%f %f (%f, %f) (%f, %f)\n",
           median,
           mean,
           etlo,
           ethi,
           hpdlo,
           hpdhi);
    #endif
    
    if (node->label)
      free(node->label);
    
    if (opt_map_hpdci)
      xasprintf(&(node->label), "[&95%%={%f, %f}]", hpdlo,hpdhi);
    else
      xasprintf(&(node->label), "[&95%%={%f, %f}]", etlo,ethi);
    node->mark = 1;

    if (opt_map_median)
      node->age = median;
    else
      node->age = mean;

    free(label);
  }

  /* set branch lengths according to ages */
  setbranchlengths(t);

  char * newick = stree_export_newick(t->root, NULL);
  fprintf(fp_out,"%s\n", newick);
  free(newick);

  stree_destroy(t,NULL);

  fclose(fp_table);
  if (opt_output)
    fclose(fp_out);
}

