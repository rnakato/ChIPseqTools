#include <stdio.h>
#include <stdlib.h>
#include "readgene.h"

int unison, nobs, maxposi, sp_scer;

const char Usage[] =
"compare_bs -1 <bs file1> -2 <bs file2>\n\
       -l: extend length (default:0) \n\
       -and: output overlapped sites \n\
       -not(default): output unique sites \n\
       -red: output overlapped sites allowing redundancy \n\
       -maxposi: use maxposition \n\
       -nobs: output only numbers \n\n";
static void argv_init(int argc, char **argv, struct peakset *, int *, int);

static int comp(const void *c1, const void *c2) {
  struct bs n1 = *(struct bs *)c1;
  struct bs n2 = *(struct bs *)c2;
  if (n1.maxIP > n2.maxIP) return -1;
  else if (n1.maxIP == n2.maxIP) return 0;
  else return 1;
}

int main(int argc, char *argv[]) {
  int i,j;
  int samplenum=2, extend_length=0;
  struct peakset *peakset = (struct peakset *)my_calloc(samplenum, sizeof(struct peakset), "peakset");
  argv_init(argc, argv, peakset, &extend_length, samplenum);

  for (i=0; i<samplenum; i++) {
    peakset[i].num = read_bs(peakset[i].file, &(peakset[i].bsarray), &(peakset[i].bsarraynum), sp_scer);
    if (peakset[i].num_defined) {
      qsort(peakset[i].bsarray, peakset[i].num, sizeof(struct bs), comp);
      peakset[i].num = min(peakset[i].num_defined,peakset[i].num);
    }
    for (j=0; j<peakset[i].num; j++) peakset[i].peakwid_total += peakset[i].bsarray[j].end - peakset[i].bsarray[j].start;
    for (j=0; j<peakset[i].num; j++) peakset[i].bsarray[j].overlap = (int *)my_calloc(samplenum, sizeof(int), "bsarray.overlap");
  }

  for (i=0; i<peakset[0].num; i++) {
    if (maxposi) {
      peakset[0].bsarray[i].start = peakset[0].bsarray[i].maxposi;
      peakset[0].bsarray[i].end   = peakset[0].bsarray[i].maxposi;
    }
  }
  for (i=0; i<samplenum; i++) {
    for (j=0; j<samplenum; j++) {
      peakset[i].bsarray_overlap[j] = (struct bs_overlap *)my_calloc(peakset[i].num+peakset[j].num, sizeof(struct bs_overlap), "bsarray_overlap");
    }
  }
  compare(peakset, 0, 1, extend_length);

  printf("#file1: %s\n#file2: %s\n", peakset[0].file, peakset[1].file);
  printf("#num1: %d\tnum2: %d\tnum1_overlap: %d (%.1f%%)\tnum1_notoverlap: %d (%.1f%%)\tnum2_overlap: %d (%.1f%%)\tnum2_notoverlap: %d (%.1f%%)\n",
         peakset[0].num, peakset[1].num,
         peakset[0].cnt_overlap[1],
         (100 * peakset[0].cnt_overlap[1] / (double)peakset[0].num),
         peakset[0].cnt_notoverlap[1],
         (100 * peakset[0].cnt_notoverlap[1] / (double)peakset[0].num),
         peakset[1].cnt_overlap[0],
         (100 * peakset[1].cnt_overlap[0] / (double)peakset[1].num),
         peakset[1].cnt_notoverlap[0],
         (100 * peakset[1].cnt_notoverlap[0] / (double)peakset[1].num));

  printf("#peakwidth total1: %d bp\tpeakwidth total2: %ul bp\toverlappeaks total: %d bp (%.2f%% / %.2f%%)\n",
         peakset[0].peakwid_total,
         peakset[1].peakwid_total,
         peakset[0].base_overlap[1],
         (peakset[0].base_overlap[1] / (double)peakset[0].peakwid_total)*100,
         (peakset[0].base_overlap[1] / (double)peakset[1].peakwid_total)*100);

  if (nobs) exit(0);

  if (unison==2) {
//    printf("chromosome\tstart\tend\tsummit 1\tsummit 2\tenrich 1\tenrich 2\tintensity 1\tintensity 2\tA (log2(1*2)/2)\tM (log2(1/2))\tdiff of summit\n");
    for (i=0; i<peakset[0].cnt_overlap_red[1]; i++) {
//      print1bs_overlap(peakset[0].bsarray_overlap[1], peakset[0].bsarray, peakset[1].bsarray, i);
      print1bsstr_overlap(peakset[0].bsarray_overlap[1], peakset[0].bsarray, peakset[1].bsarray, i);
    }
  }else if (unison==1) {
//   printf("chromosome\tstart\tend\tsummit\tlength\tenrich\tintensity\n");
    for (i=0; i<peakset[0].num; i++) {
      if (peakset[0].bsarray[i].overlap[1]) {
//        print1bs(peakset[0].bsarray, i);
        print1bsstr(peakset[0].bsarray, i);
      }
    }
  }else{
//    printf("chromosome\tstart\tend\tsummit\tlength\tenrich\tintensity\n");
    for (i=0; i<peakset[0].num; i++) {
      if (!(peakset[0].bsarray[i].overlap[1])) {
 //       print1bs(peakset[0].bsarray, i);
        print1bsstr(peakset[0].bsarray, i);
      }
    }
  }


  for (i=0; i<samplenum; i++) {
    free(peakset[i].bsarray);
    for (j=0; j<samplenum; j++) free(peakset[i].bsarray_overlap[j]);
    free(peakset[i].bsarray_overlap);
    free(peakset[i].base_overlap);
    free(peakset[i].cnt_overlap);
    free(peakset[i].cnt_notoverlap);
    free(peakset[i].cnt_overlap_red);
  }
  return 0;
}

static void argv_init(int argc, char **argv, struct peakset *peakset, int *extend_length, int samplenum)
{
  int i;
  unison=0;
  maxposi=0;
  nobs=0;
  sp_scer=0;
  for (i=0; i<samplenum; i++) {
    peakset[i].file = NULL;
    peakset[i].num = 0;
    peakset[i].num_defined = 0;
    peakset[i].bsarraynum = STRUCT_BS_MAX;
    peakset[i].bsarray = (struct bs *)my_calloc(peakset[i].bsarraynum, sizeof(struct bs), "bsarray");
    peakset[i].bsarray_overlap = (struct bs_overlap **)my_calloc(samplenum, sizeof(struct bs_overlap *), "*bsarray_overlap");
    peakset[i].peakwid_total = 0;
    peakset[i].cnt_overlap = (int *)my_calloc(samplenum, sizeof(int), "cnt_overlap");
    peakset[i].cnt_notoverlap = (int *)my_calloc(samplenum, sizeof(int), "cnt_notoverlap");
    peakset[i].cnt_overlap_red = (int *)my_calloc(samplenum, sizeof(int), "cnt_overlap_red");
    peakset[i].base_overlap = (int *)my_calloc(samplenum, sizeof(int), "base_overlap");
  }

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i], "-1"))              peakset[0].file = argv[++i];
    else if (!strcmp(argv[i], "-2"))         peakset[1].file = argv[++i];
    else if (!strcmp(argv[i], "-n1"))        peakset[0].num_defined = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-n2"))        peakset[1].num_defined = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-l"))         *extend_length=atoi(argv[++i]);
    else if (!strcmp(argv[i], "-red"))       unison=2;
    else if (!strcmp(argv[i], "-and"))       unison=1;
    else if (!strcmp(argv[i], "-not"))       unison=0;
    else if (!strcmp(argv[i], "-maxposi"))   maxposi=1;
    else if (!strcmp(argv[i], "-scer"))      sp_scer=1;
    else if (!strcmp(argv[i], "-nobs"))      nobs=1;
    else goto err;
  }
  if (!(peakset[0].file) || !(peakset[1].file)) goto err;
  return;

 err:
  fprintf(stderr, "%s", Usage);
  exit(0);
}
