#include <stdio.h>
#include <stdlib.h>
#include "readgene.h"

static void cnt_overlap(struct peakset *peakset, int sample1, int sample2)
{
  int i;
  for (i=0; i<peakset[sample1].num; i++){
    if (peakset[sample1].bsarray[i].overlap[sample2]){
      peakset[sample1].cnt_overlap[sample2]++;
    }else{
      peakset[sample1].cnt_notoverlap[sample2]++;
    }
  }
}

static void copy_bs2bs_overlap(struct peakset *peakset, int sample1, int sample2, int i, int j, int num)
{
  peakset[sample1].bsarray_overlap[sample2][num].bs1 = i;
  peakset[sample1].bsarray_overlap[sample2][num].bs2 = j;
  if (peakset[sample2].bsarray[j].enrich) peakset[sample1].bsarray_overlap[sample2][num].ratio_enrich = peakset[sample1].bsarray[i].enrich/peakset[sample2].bsarray[j].enrich;
  else peakset[sample1].bsarray_overlap[sample2][num].ratio_enrich = 0;
  if (peakset[sample2].bsarray[j].maxIP) peakset[sample1].bsarray_overlap[sample2][num].ratio_maxIP = peakset[sample1].bsarray[i].maxIP/peakset[sample2].bsarray[j].maxIP;
  else peakset[sample1].bsarray_overlap[sample2][num].ratio_maxIP = 0;
  peakset[sample1].bsarray_overlap[sample2][num].peaksummitdiff = abs(peakset[sample1].bsarray[i].maxposi - peakset[sample2].bsarray[j].maxposi);
}

static void getoverlap(struct peakset *peakset, int sample1, int sample2, int i, int j, int num)
{
  peakset[sample1].bsarray[i].overlap[sample2]=1;
  copy_bs2bs_overlap(peakset, sample1, sample2, i, j, num);
  peakset[sample1].base_overlap[sample2] += min(peakset[sample1].bsarray[i].end, peakset[sample2].bsarray[j].end) - max(peakset[sample1].bsarray[i].start, peakset[sample2].bsarray[j].start);
}

void compare(struct peakset *peakset, int sample1, int sample2, int extend_length)
{
  int i,j, num=0;
  for (i=0; i<peakset[sample1].num; i++){
    for (j=0; j<peakset[sample2].num; j++){
      if (!strcmp(peakset[sample1].bsarray[i].chr, peakset[sample2].bsarray[j].chr)
         && peakset[sample2].bsarray[j].start <= peakset[sample1].bsarray[i].end + extend_length
         && peakset[sample2].bsarray[j].end >= peakset[sample1].bsarray[i].start - extend_length){
        getoverlap(peakset, sample1, sample2, i, j, num);
        if (sample1 != sample2) getoverlap(peakset, sample2, sample1, j, i, num);
        num++;
      }
    }
  }
  cnt_overlap(peakset, sample1, sample2);
  if (sample1 != sample2) cnt_overlap(peakset, sample2, sample1);
  peakset[sample1].cnt_overlap_red[sample2] = num;
  return;
}
