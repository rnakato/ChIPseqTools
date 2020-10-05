#include <stdio.h>
#include <stdlib.h>
#include "readgene.h"
#define TOPNUM 3000
#define NUM_CRM 100000
#define SAMPLE_NUM_MAX 1024

int n_toprank, sp_scer, quantitative;

const char Usage[] =
"merge_bs2CRM -i <bs file> -name <name> [-i <bs file> -name <name> ...]\n\
       -l: extend length (default:0) \n\
       -n: number of peaks for clustering (default:3000, setting 0 means use all peaks) \n\
       -qnt: quantitative analysis \n\n";

static int argv_init(int argc, char **argv, struct peakset **refpeakset, int *extend_length);
static struct peakset *peakset_new(int IPnum, char *file[], char *name[]);

struct CRM{
  char chr[128];
  int start;
  int end;
  double occur[SAMPLE_NUM_MAX];
};

static int comp(const void *c1, const void *c2) {
  struct bs n1 = *(struct bs *)c1;
  struct bs n2 = *(struct bs *)c2;
  if      (n1.maxIP > n2.maxIP)  return -1;
  else if (n1.maxIP == n2.maxIP) return 0;
  else return 1;
}

static void merge_overlappedCRM(struct CRM *CRM, int samplenum, int *num_CRM) {
  int i,j,k;
  for (i=0; i<*num_CRM; i++) {
    for (j=i+1; j<*num_CRM; j++) {
      if (!strcmp(CRM[j].chr, CRM[i].chr)
	  && CRM[j].start <= CRM[i].end
	  && CRM[j].end   >= CRM[i].start)
	{
	  // jをiにマージ, jにはラストのオブジェクトを代入
	  CRM[i].start = min(CRM[i].start, CRM[j].start);
	  CRM[i].end   = max(CRM[i].end,   CRM[j].end);
	  for (k=0; k<samplenum; k++) {
	    CRM[i].occur[k] = max(CRM[i].occur[k], CRM[j].occur[k]);
	  }
	  CRM[j] = CRM[*num_CRM -1];
	  (*num_CRM)--;
	}
    }
  }
  return;
}


static void renew_CRM(struct CRM *CRM, int k, struct bs *bsarray, int i, int j) {
  CRM[k].start = min(CRM[k].start, bsarray[j].start);
  CRM[k].end   = max(CRM[k].end,   bsarray[j].end);
  if (quantitative) CRM[k].occur[i] = bsarray[j].maxIP; // peak intensity
  else              CRM[k].occur[i] = 1;                // binary 0 or 1
}

static void add_CRM(struct CRM *CRM, int k, struct bs *bsarray, int i, int j) {
  strcpy(CRM[k].chr, bsarray[j].chr);
  CRM[k].start = bsarray[j].start;
  CRM[k].end   = bsarray[j].end;
  if (quantitative) CRM[k].occur[i] = bsarray[j].maxIP; // peak intensity
  else              CRM[k].occur[i] = 1;                // binary 0 or 1
}

int main(int argc, char *argv[]) {
  int i,j,k, extend_length=0;
  struct peakset *peakset = NULL;
  int samplenum = argv_init(argc, argv, &peakset, &extend_length);
  for (i=0; i<samplenum; i++) {
    peakset[i].num = read_bs(peakset[i].file, &(peakset[i].bsarray), &(peakset[i].bsarraynum), sp_scer);
    if (n_toprank && n_toprank < peakset[i].num) {
      qsort(peakset[i].bsarray, peakset[i].num, sizeof(struct bs), comp);
      peakset[i].num = n_toprank;
    }
  }

  int num_CRM = 0;
  int num_CRMarray = NUM_CRM;
  struct CRM *CRM = (struct CRM *)my_calloc(num_CRMarray, sizeof(struct CRM), "CRM");
  int on;

  for (i=0; i<samplenum; i++) {
    for (j=0; j<peakset[i].num; j++) {
      on=0;
      for (k=0; k<num_CRM; k++) {
	if (!strcmp(peakset[i].bsarray[j].chr, CRM[k].chr)
	   && peakset[i].bsarray[j].start <= CRM[k].end + extend_length
	   && peakset[i].bsarray[j].end >= CRM[k].start - extend_length) {
	  renew_CRM(CRM, k, peakset[i].bsarray, i, j);
	  on=1;
	}
      }
      if (!on) {
	add_CRM(CRM, num_CRM, peakset[i].bsarray, i, j);
	num_CRM++;
	if (num_CRM >= num_CRMarray) {
	  num_CRMarray += NUM_CRM;
	  CRM = (struct CRM *)my_realloc(CRM, num_CRMarray*sizeof(struct CRM), "CRM");
	}
      }
    }
    // 得られた（伸長された）CRMで重なるものをマージ
    merge_overlappedCRM(CRM, samplenum, &num_CRM);
  }

  printf("chromosome\tstart\tend\tlength\tID\tnum of peaks");
  for (i=0; i<samplenum; i++) {
    printf("\t%s", peakset[i].name);
    if (i==samplenum-1) printf("\n");
  }
  int num;
  for (i=0; i<num_CRM; i++) {
    num=0;
    for (j=0; j<samplenum; j++) {
      if (CRM[i].occur[j]>0) num++;
    }
    printf("chr%s\t%d\t%d\t%d\tCRM%d\t%d", CRM[i].chr, CRM[i].start, CRM[i].end, CRM[i].end - CRM[i].start, i, num);
    for (j=0; j<samplenum; j++) {
      if (quantitative) printf("\t%f",   CRM[i].occur[j]);
      else              printf("\t%.0f", CRM[i].occur[j]);
      if (j==samplenum-1) printf("\n");
    }
  }

  for (i=0; i<samplenum; i++) {
    free(peakset[i].bsarray);
    for (j=0; j<samplenum; j++) free(peakset[i].bsarray_overlap[j]);
    free(peakset[i].bsarray_overlap);
    free(peakset[i].base_overlap);
    free(peakset[i].cnt_overlap);
    free(peakset[i].cnt_notoverlap);
  }
  return 0;
}

static int argv_init(int argc, char **argv, struct peakset **refpeakset, int *extend_length) {
  int i, IPnum=0, namenum=0;
  char *file[100], *name[100];
  sp_scer=0;
  n_toprank = TOPNUM;
  quantitative = 0;

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i], "-i"))              file[IPnum++] = argv[++i];
    else if (!strcmp(argv[i], "-qnt"))       quantitative = 1;
    else if (!strcmp(argv[i], "-name"))      name[namenum++] = argv[++i];
    else if (!strcmp(argv[i], "-l"))         *extend_length = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-n"))         n_toprank = atoi(argv[++i]);
    else {
      fprintf(stderr, "Invalid argument: %s\n", argv[i]);
      goto err;
    }
  }
  if (!IPnum) {
    fprintf(stderr, "no input file.\n");
    goto err;
  }
  if (IPnum != namenum) {
    fprintf(stderr, "inconsistent number: -i and -name.\n");
    goto err;
  }
  if (IPnum >= SAMPLE_NUM_MAX) {
    fprintf(stderr, "error: input file num is >%d.\n", SAMPLE_NUM_MAX);
    goto err;
  }

  *refpeakset = peakset_new(IPnum, file, name);
  return IPnum;

 err:
  fprintf(stderr, "%s", Usage);
  exit(0);
}

static struct peakset *peakset_new(int IPnum, char *file[], char *name[]) {
  int i;
  struct peakset *peakset = (struct peakset *)my_calloc(IPnum, sizeof(struct peakset), "peakset");
  for (i=0; i<IPnum; i++) {
    peakset[i].file = file[i];
    peakset[i].name = name[i];
    peakset[i].num = 0;
    peakset[i].bsarraynum = STRUCT_BS_MAX;
    peakset[i].bsarray = (struct bs *)my_calloc(peakset[i].bsarraynum, sizeof(struct bs), "bsarray");
    peakset[i].bsarray_overlap = (struct bs_overlap **)my_calloc(IPnum, sizeof(struct bs_overlap *), "*bsarray_overlap");
    peakset[i].peakwid_total = 0;
    peakset[i].cnt_overlap = (int *)my_calloc(IPnum, sizeof(int), "cnt_overlap");
    peakset[i].cnt_notoverlap = (int *)my_calloc(IPnum, sizeof(int), "cnt_notoverlap");
    peakset[i].base_overlap = (int *)my_calloc(IPnum, sizeof(int), "base_overlap");
  }
  return peakset;
}
