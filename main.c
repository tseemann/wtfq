#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <getopt.h>
#include <math.h>
#include <assert.h>

//#define DEBUG

// Heng's FASTQ reader
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// software attributes
#define EXENAME "wtfq"
#define GITHUB_URL "https://github.com/tseemann/wtfq"
#define AUTHOR "Torsten Seemann"
#define VERSION "0.1"

// default parameter values
#define DEFAULT_CORRECTION 1.5
#define DEFAULT_SEED 42

// warning: implicit declaration of function 'fileno' [-Wimplicit-function-declaration]
int fileno(FILE *);
void srand48(long int);
double drand48(void);

//------------------------------------------------------------------------
void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  fprintf(out, "SYNOPSIS\n  Generate extra FASTQ reads for extreme-GC/undersequenced regions\n");
  fprintf(out, "USAGE\n");
  fprintf(out, "  %s [options] reads.fq.gz > extras.fq\n", EXENAME);
  fprintf(out, "  seqtk mergepe R1.fq.gz R2.fq.gz | %s - | gzip > extras.fq.gz \n", EXENAME);
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -h        Show this help\n");
  fprintf(out, "  -v        Print version and exit\n");
  fprintf(out, "  -q        Quiet mode; not non-error output\n");
  fprintf(out, "  -c %5.3f  Correction exponent\n", DEFAULT_CORRECTION);
  fprintf(out, "  -s %-4d   Random seed\n", DEFAULT_SEED);
//  fprintf(out, "  -p STRING  Prefix if not\n");
//  fprintf(out, "  -o FILE1   Output file for R1\n");
//  fprintf(out, "  -O FILE2   Output file for R2\n");
  fprintf(out, "URL\n  %s (%s)\n", GITHUB_URL, AUTHOR);
  exit(retcode);
}

//------------------------------------------------------------------------
double kseq_gc(const kseq_t* r) {
  char* s = r->seq.s;
  size_t GC=0, N=0, Nbad=0;
  
  while (*s++) {
    switch (toupper(*s)) {
      case 'A' : 
      case 'T' : 
        N++; break;
      case 'G' : 
      case 'C' : 
        N++; GC++; break;
      default : 
        Nbad++; // ambiguous nucleotides
    }
  }

#ifdef DEBUG
    fprintf(stderr, "# %s\n# stats: GC=%lu/%lu ignored=%lu/%lu\n\n", 
      r->seq.s, GC, N, Nbad, r->seq.l);
#endif
  
  return N==0 ? -1.0 : (GC / (double) N);
}

//------------------------------------------------------------------------

void kseq_write_dupes(const kseq_t* s, int num, FILE* stream) 
{
  static int counter=0;
  for (int i=0; i < num; i++) {
    // if the length of the quality string is zero, it was fasta input
    fputc( s->qual.l ? '@' : '>', stream);
    fprintf(stream, "%s:%s-%d\n%s\n", s->name.s, EXENAME, counter++, s->seq.s);
    if (s->qual.l) fprintf(stream, "+\n%s\n", s->qual.s);
  }
}
 
//------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // parse command line parameters
  int opt, quiet=0;
  double corr=DEFAULT_CORRECTION;
  long seed=DEFAULT_SEED;
  
  while ((opt = getopt(argc, argv, "hvqfc:s:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet=1; break;
      case 'c': corr=atof(argv[optind]); break;
      case 's': seed=atol(argv[optind]); break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      default : show_help(EXIT_FAILURE);
    }
  } 

  // correction factor must be in [0,1]
  if (corr < 1.0) {
    fprintf(stderr, "ERROR: bad -c %f ; must be in [0,Inf)\n", corr);
    exit(EXIT_FAILURE);
  }

  // random seed
  srand48(seed);
  if (!quiet) fprintf(stderr, "RUNNING : %s %s\n", EXENAME, VERSION);

  // open stdin or the file on commandline
  FILE* input = NULL;
  if (optind >= argc) {
    // no file specified so check if stdin is a pipe
    if (isatty(fileno(stdin))) {
      show_help(EXIT_FAILURE);
    }
    else {
      input = stdin;
    }
  }
  else {
    // opening a filename
    const char* seqfn = argv[optind];
    // "-" is the special filename for stdin
    input = strcmp(seqfn, "-") ? fopen(seqfn, "r") : stdin;
    if (! input) {
      fprintf(stderr, "ERROR: Could not open '%s'\n", seqfn);
      exit(EXIT_FAILURE);
    }
    if (!quiet) fprintf(stderr, "READING : %s\n", input==stdin ? "stdin" : seqfn);
  }

  // open filehandle with zlib
  gzFile fp = gzdopen(fileno(input), "r"); 
  if (! fp) {
    fprintf(stderr, "ERROR: Could not gzopen input\n");
    exit(EXIT_FAILURE);  
  }
  
  long l, N=0, L=0, D=0;
  double gc_sum=0;
  kseq_t* kseq = kseq_init(fp);
 
  while ((l = kseq_read(kseq)) >= 0) {    

    int copies=0;
    double gc = kseq_gc(kseq);
    if (gc >= 0) {  // no AGTCs in read   
      gc_sum += gc;
      double GC = gc_sum / (double) (N+1) ;
      double dev = -(gc-GC)/GC;
      if (dev > 0) {
        copies = (int) (pow(GC/gc,corr) + 0.5);  // invented formula w/ no evidence
        if (copies <= 1) copies=0;
      }
    }
 
    if (copies > 0) {
      kseq_write_dupes(kseq, copies, stdout);
      D += copies;
    }

    // keep stats
    L += kseq->seq.l;
    N++;
  }
  kseq_destroy(kseq); 
  gzclose(fp); 

  // reveal status
  if (!quiet) {
    fprintf(stderr, "INPUT   : seqs=%ld bp=%ld avglen=%ld gc=%lf\n", N, L, L/N, gc_sum/N);
    fprintf(stderr, "OUTPUT  : seqs=%ld\n", D);
  }

  return 0;
}

//------------------------------------------------------------------------

