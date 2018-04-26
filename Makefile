
EXE = wtfq
PREFIX = /usr/local
TESTDIR = test

CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c99 
LIBS = -lz -lm

.PHONY: check clean install
.DEFAULT: all

all: $(EXE)

$(EXE): main.c 
	$(CC) $(CFLAGS) -o $(EXE) $^ $(LIBS)

main.c: kseq.h

install: $(EXE)
	install -v -t $(PREFIX)/bin $(EXE)

clean:
	$(RM) *~ *.o $(EXE)

test:
	./$(EXE) -v
	./$(EXE) -h
#	./$(EXE) $(TESTDIR)/R1.fq.gz $(TESTDIR)/R2.fq.gz

