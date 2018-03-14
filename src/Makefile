# Copyright (C) 2018 Tomas Flouri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
# Department of Genetics, Evolution and Environment,
# University College London, Gower Street, London WC1E 6BT, United Kingdom

# Profiling
#PROFILING=-g -pg
PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare

ifndef CC
CC = gcc
endif
CFLAGS = -g -O3 -D_GNU_SOURCE $(WARN)
LINKFLAGS=$(PROFILING)
LIBS=-lm

PROG=summarizer

all: $(PROG)

OBJS=summarizer.o summary.o util.o arch.o combine.o parse.o summaryfull.o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $+ $(LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
