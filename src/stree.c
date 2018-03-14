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

static char * stree_export_newick_recursive(const snode_t * root,
                                            char * (*cb_serialize)(const snode_t *))
{
  char * newick;
  int size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = stree_export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
    {
      return NULL;
    }
    char * subtree2 = stree_export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s",
                               subtree1,
                               subtree2,
                               temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s:%f",
                               subtree1,
                               subtree2,
                               root->label ? root->label : "",
                               root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
  {
    fatal("Memory allocation during newick export failed");
  }

  return newick;
}

char * stree_export_newick(const snode_t * root,
                           char * (*cb_serialize)(const snode_t *))
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = stree_export_newick_recursive(root->left,cb_serialize);
    if (!subtree1)
      fatal("Unable to allocate enough memory.");

    char * subtree2 = stree_export_newick_recursive(root->right,cb_serialize);
    if (!subtree2)
      fatal("Unable to allocate enough memory.");

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s;",
                               subtree1,
                               subtree2,
                               temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s:%f;",
                               subtree1,
                               subtree2,
                               root->label ? root->label : "",
                               root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("memory allocation during newick export failed");

  return newick;
}

