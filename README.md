# MCMCtree summarizer

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Compilation instructions

The program can be compiled using the command:

```bash
make
```

## Running summarizer

You can run the summarizer using the command:

```bash
summarizer --summarize FILENAME --output OUTFILE
```

where `FILENAME` is an MCMCtree MCMC sample file (typically `mcmc.txt`).
Summary is written in OUTFILE.

You can also combine several MCMC files using the command:

```bash
summarizer --combine FILENAME --output OUTFILE
```

where `FILENAME` is a file containing a list of files to be combined. The
output is written in OUTFILE. The program creates one additional file called
OUTFILE.index, and fills it with the number of samples in each of the files to
be combined. This file can be later passed to the `--summarize` command in the
following way (as INDEXFILE):

```bash
summarizer --summarize FILENAME --output OUTFILE --index INDEXFILE
```

in order to create per-dataset summaries along with the combined summary. The
individual dataset summaries are written in separate files as OUTFILE.1.txt,
OUTFILE.table.1.txt, OUTFILE.2.txt, OUTFILE.table.2.txt, ..., OUTFILE.N.txt
OUTFILE.table.N.txt, where **N** is the number of combined datasets. The combined
summary is written as OUTFILE.combined.txt and OUTFILE.table.combined.txt.

Please note that in both cases (`--summarize and --combine)`, summarixer
ignores the first line of each processed file, which is typically the header
line containing the column labels. If you wish to ignore more lines, please use
the option `--skip INTEGER` (default: 1).

Finally, summarizer can map the table summaries to a target tree. The command is:

```bash
summarizer --map FILENAME --tree TREEFILE --output OUTFILE
```

This will read the newick tree from file TREEFILE (internal nodes must have numeric labels) and will 
create branch lengths equal to the posterior mean age difference between nodes. Furthermore, it will
add the Equal-Tail CI as attribute to inner nodes, such that the resulting tree can be viewed using
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/). The tree is saved in OUTFILE (if specified)
or printed on screen otherwise.
If you wish to use the median age instead of mean ages, supply the argument `--median`. If you wish
to use the 95% HPD CI instead of Equal-tail CI, supply the argument `--hpdci`.
