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
output is written in OUTFILE.

Please note that in both cases (`--summarize and --combine)`, summarixer
ignores the first line of each processed file, which is typically the header
line containing the column labels. If you wish to ignore more lines, please use
the option `--skip INTEGER` (default: 1).

```
