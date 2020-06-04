#ifndef _READDATA_H_
#define _READDATA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "my.h"
#include "stringp.h"

#define STRUCT_BS_MAX 100000
#define STRUCT_GENE_MAX 100000

typedef enum{CODINGGENE, PSEUDO, RNA, OTHERS, PROCESS} genetype;
typedef enum{NONE, UPSTREAM, DOWNSTREAM} updown;

typedef enum{
  REFSEQ,
  ENSEMBL,
  SCER_GENE,
  POMBE_GENE,
  SCHYZON_GENE
} genefiletype;

struct bs{
  char line[10240];
  char chr[128];
  int start;
  int end;
  int maxposi;
  double enrich;
  double maxIP;
  int *overlap;
  int totss;
  updown updown_fromtss;
  int geneid;
};

struct bs_overlap{
  int bs1, bs2;
  double ratio_enrich;
  double ratio_maxIP;
  int peaksummitdiff;
};

struct peakset{
  char *file;
  char *name;
  int num;
  int num_defined;
  struct bs *bsarray;
  struct bs_overlap **bsarray_overlap;
  int bsarraynum;
  int peakwid_total;
  int *base_overlap;
  int *cnt_overlap;
  int *cnt_notoverlap;
  int *cnt_overlap_red;
};

struct exon{
  int start;
  int end;
};

struct gene{
  char name[128];
  char type[32];
  char desc[1028];
  char ID[32];
  int start;
  int end;
  int exonnum;
  struct exon *exon;
  char chr[128];
  int dir;
  genetype genetype;
  int genic;
  int up;
  int down;
  int on;
  int bsid_up;
  int bsid_down;
};

struct genename{
  char name[128];
  int dist;
};

int read_bs(char *argv, struct bs **bsarray, int *bsarraynum, int);
int read_genelist(char *genefile, struct genename **genelist);
int read_gene_ENS(char *filename, struct gene **ref_gene, int);
int read_refFlat(char *filename, struct gene **ref_gene, int);
int read_refseq(char *filename, struct gene **ref_gene);
void print1bs(struct bs *bsarray, int i);
void print1bsstr(struct bs *bsarray, int i);
void print1bs_overlap(struct bs_overlap *bsarray, struct bs *, struct bs *, int i);
void print1bsstr_overlap(struct bs_overlap *bsarray, struct bs *, struct bs *, int i);

void compare(struct peakset *peakset, int sample1, int sample2, int extend_length);

#define max(a, b) (((a) > (b))?(a) :(b))
#define min(a, b) (((a) < (b))?(a) :(b))

#endif
