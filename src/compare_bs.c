#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int unison, nobs, maxposi, sp_scer, showallcols;

#define MYFREE(p) {if(p){free(p); (p)=NULL;} }
#define max(a, b) (((a) > (b))?(a) :(b))
#define min(a, b) (((a) < (b))?(a) :(b))

#define STR_LEN 100000
#define ELEM_NUM 256
#define ELEM_NUM2 2560

struct elem{
  char str[8000];
};

#define STRUCT_BS_MAX 100000
typedef enum{NONE, UPSTREAM, DOWNSTREAM} updown;

struct bs{
  char line[12800];
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

const char Usage[] =
"compare_bs -1 <bs file1> -2 <bs file2>\n\
       -l: extend length (default:0) \n\
       -and: output overlapped sites \n\
       -not(default): output unique sites \n\
       -red: output overlapped sites allowing redundancy \n\
       -maxposi: use maxposition \n\
       -showallcols: show all columns of the first input file \n\
       -nobs: output only numbers \n\n";
static void argv_init(int argc, char **argv, struct peakset *, int *, int);

static int comp(const void *c1, const void *c2){
  struct bs n1 = *(struct bs *)c1;
  struct bs n2 = *(struct bs *)c2;
  if(n1.maxIP > n2.maxIP) return -1;
  else if(n1.maxIP == n2.maxIP) return 0;
  else return 1;
}

FILE *my_fopen_r(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "r"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

void *my_calloc(size_t n, size_t s, char *name){
  void *p;
  p = calloc(n,s);
  if(!p){
    fprintf(stderr,"[E]failed calloc: %s\n", name); 
    exit(1);
  }
  return p;
}

void *my_realloc(void *p, size_t s, char *name){
  p = realloc(p,s);
  if(!p){
    fprintf(stderr,"[E]failed realloc: %s\n", name); 
    exit(1);
  }
  return p;
}

void chomp(char *str){
  char *p = strchr(str, '\n');
  if(p) p[0]='\0';
  return;
}

int ParseLine(char *str, struct elem clm[]){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len, sizeof(char), "ParseLine");
  for(i=0; i<=len; i++){
    if(str[i]=='\0' || str[i]=='\n'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      MYFREE(strtemp);
      return ++num;
    }
    if(str[i]=='\t'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++; 
      if(num >= ELEM_NUM){
	fprintf(stderr, "error: too many columns: %s", str);
	exit(0);
      }
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  MYFREE(strtemp);
  return num;
}

static char *checkchrname(char *name){
  char *p = strstr(name, "chr");
  if(p) return p+3;
  else return name;
}


static void changechr_yeast(struct bs **bsarray, int num, char *str, int sp_scer){
  char *chrstr= checkchrname(str);
  if(!strcmp(chrstr, "I"))         strcpy((*bsarray)[num].chr, "1");
  else if(!strcmp(chrstr, "II"))   strcpy((*bsarray)[num].chr, "2");
  else if(!strcmp(chrstr, "III"))  strcpy((*bsarray)[num].chr, "3");
  else if(!strcmp(chrstr, "IV"))   strcpy((*bsarray)[num].chr, "4");
  else if(!strcmp(chrstr, "V"))    strcpy((*bsarray)[num].chr, "5");
  else if(!strcmp(chrstr, "VI"))   strcpy((*bsarray)[num].chr, "6");
  else if(!strcmp(chrstr, "VII"))  strcpy((*bsarray)[num].chr, "7");
  else if(!strcmp(chrstr, "VIII")) strcpy((*bsarray)[num].chr, "8");
  else if(!strcmp(chrstr, "IX"))   strcpy((*bsarray)[num].chr, "9");
  else if(!strcmp(chrstr, "X")){
    if(sp_scer) strcpy((*bsarray)[num].chr, "10"); else strcpy((*bsarray)[num].chr, "X");
  }else if(!strcmp(chrstr, "XI"))   strcpy((*bsarray)[num].chr, "11");
  else if(!strcmp(chrstr, "XII"))  strcpy((*bsarray)[num].chr, "12");
  else if(!strcmp(chrstr, "XIII")) strcpy((*bsarray)[num].chr, "13");
  else if(!strcmp(chrstr, "XIV"))  strcpy((*bsarray)[num].chr, "14");
  else if(!strcmp(chrstr, "XV"))   strcpy((*bsarray)[num].chr, "15");
  else if(!strcmp(chrstr, "XVI"))  strcpy((*bsarray)[num].chr, "16");
  else strcpy((*bsarray)[num].chr, chrstr);
  return;
}

static void copy_bs2bs_overlap(struct peakset *peakset, int sample1, int sample2, int i, int j, int num){
  peakset[sample1].bsarray_overlap[sample2][num].bs1 = i;
  peakset[sample1].bsarray_overlap[sample2][num].bs2 = j;
  if(peakset[sample2].bsarray[j].enrich) peakset[sample1].bsarray_overlap[sample2][num].ratio_enrich = peakset[sample1].bsarray[i].enrich/peakset[sample2].bsarray[j].enrich;
  else peakset[sample1].bsarray_overlap[sample2][num].ratio_enrich = 0;
  if(peakset[sample2].bsarray[j].maxIP) peakset[sample1].bsarray_overlap[sample2][num].ratio_maxIP = peakset[sample1].bsarray[i].maxIP/peakset[sample2].bsarray[j].maxIP;
  else peakset[sample1].bsarray_overlap[sample2][num].ratio_maxIP = 0;
  peakset[sample1].bsarray_overlap[sample2][num].peaksummitdiff = abs(peakset[sample1].bsarray[i].maxposi - peakset[sample2].bsarray[j].maxposi);
}

int read_bs(char *filename, struct bs **bsarray, int *bsarraynum, int sp_scer){
  int num=0, linenum;
  char *chrstr=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  struct elem clm[ELEM_NUM];

  FILE *IN = my_fopen_r(filename);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n' || str[0]=='#') continue;
    chomp(str);

    linenum = ParseLine(str, clm);
    if(linenum<3) continue;
    if(strstr(clm[0].str, "num1")) continue;
    if(!strcmp(clm[0].str, "chromosome")) continue;

    chrstr = checkchrname(clm[0].str);
    changechr_yeast(bsarray, num, chrstr, sp_scer);
    (*bsarray)[num].start = atoi(clm[1].str);
    (*bsarray)[num].end = atoi(clm[2].str);
    if(clm[3].str){
      (*bsarray)[num].maxposi = atoi(clm[3].str);
      if(atoi(clm[3].str)< (*bsarray)[num].start || atoi(clm[3].str)> (*bsarray)[num].end) (*bsarray)[num].maxposi = ((*bsarray)[num].start + (*bsarray)[num].end)/2;
    }
    if(linenum>=6) (*bsarray)[num].enrich  = atof(clm[5].str);
    if(linenum>=7) (*bsarray)[num].maxIP   = atof(clm[6].str);
    strcpy((*bsarray)[num].line, str);
    num++;
    if(num>=*bsarraynum){
      *bsarraynum += STRUCT_BS_MAX;
      *bsarray = (struct bs *)my_realloc(*bsarray, *bsarraynum * sizeof(struct bs), "bsarray");
    }
  }
  fclose(IN);
  free(str);
  return num;
}

static void getoverlap(struct peakset *peakset, int sample1, int sample2, int i, int j, int num, int extend_length){
  peakset[sample1].bsarray[i].overlap[sample2]=1;
  copy_bs2bs_overlap(peakset, sample1, sample2, i, j, num);
  peakset[sample1].base_overlap[sample2] += min(peakset[sample1].bsarray[i].end, peakset[sample2].bsarray[j].end) - max(peakset[sample1].bsarray[i].start, peakset[sample2].bsarray[j].start) + extend_length;
}

static void cnt_overlap(struct peakset *peakset, int sample1, int sample2){
  int i;
  for(i=0; i<peakset[sample1].num; i++){
    if(peakset[sample1].bsarray[i].overlap[sample2]){
      peakset[sample1].cnt_overlap[sample2]++;
    }else{
      peakset[sample1].cnt_notoverlap[sample2]++;
    }
  }
}

void compare(struct peakset *peakset, int sample1, int sample2, int extend_length){
  int i,j, num=0;
  for(i=0; i<peakset[sample1].num; i++){
    for(j=0; j<peakset[sample2].num; j++){
      if(!strcmp(peakset[sample1].bsarray[i].chr, peakset[sample2].bsarray[j].chr) 
	 && peakset[sample2].bsarray[j].start <= peakset[sample1].bsarray[i].end + extend_length 
	 && peakset[sample2].bsarray[j].end >= peakset[sample1].bsarray[i].start - extend_length){
	getoverlap(peakset, sample1, sample2, i, j, num, extend_length);
	if(sample1 != sample2) getoverlap(peakset, sample2, sample1, j, i, num, extend_length);
	num++;
      }
    }
  }
  cnt_overlap(peakset, sample1, sample2);
  if(sample1 != sample2) cnt_overlap(peakset, sample2, sample1);
  peakset[sample1].cnt_overlap_red[sample2] = num;
  return;
}

void print1bs(struct bs *bsarray, int i){
  if(showallcols) printf("%s\n", bsarray[i].line);
  else printf("chr%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", bsarray[i].chr, bsarray[i].start, bsarray[i].end, bsarray[i].maxposi, bsarray[i].end-bsarray[i].start, bsarray[i].enrich, bsarray[i].maxIP);
}

void print1bs_overlap(struct bs_overlap *bsarray, struct bs *bs1, struct bs *bs2, int i){
  int s1 = bsarray[i].bs1;
  int s2 = bsarray[i].bs2;
  printf("chr%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", bs1[s1].chr, bs1[s1].start, bs1[s1].end, bs1[s1].maxposi, bs2[s2].maxposi, bs1[s1].enrich, bs2[s2].enrich, bs1[s1].maxIP, bs2[s2].maxIP, log2(bs1[s1].maxIP*bs2[s2].maxIP)/2, log2(bsarray[i].ratio_maxIP), bsarray[i].peaksummitdiff);
}


int main(int argc, char *argv[]){
  int i,j;
  int samplenum=2, extend_length=0;
  struct peakset *peakset = (struct peakset *)my_calloc(samplenum, sizeof(struct peakset), "peakset");
  argv_init(argc, argv, peakset, &extend_length, samplenum);

  for(i=0; i<samplenum; i++){
    peakset[i].num = read_bs(peakset[i].file, &(peakset[i].bsarray), &(peakset[i].bsarraynum), sp_scer);
    if(peakset[i].num_defined){
      qsort(peakset[i].bsarray, peakset[i].num, sizeof(struct bs), comp);
      peakset[i].num = min(peakset[i].num_defined,peakset[i].num);
    }
    for(j=0; j<peakset[i].num; j++) peakset[i].peakwid_total += peakset[i].bsarray[j].end - peakset[i].bsarray[j].start; 
    for(j=0; j<peakset[i].num; j++) peakset[i].bsarray[j].overlap = (int *)my_calloc(samplenum, sizeof(int), "bsarray.overlap");
  }

  for(i=0; i<peakset[0].num; i++){
    if(maxposi){
      peakset[0].bsarray[i].start = peakset[0].bsarray[i].maxposi;
      peakset[0].bsarray[i].end   = peakset[0].bsarray[i].maxposi;   
    }
  }
  for(i=0; i<samplenum; i++){
    for(j=0; j<samplenum; j++){
      peakset[i].bsarray_overlap[j] = (struct bs_overlap *)my_calloc(peakset[i].num+peakset[j].num, sizeof(struct bs_overlap), "bsarray_overlap");
    }
  }
  compare(peakset, 0, 1, extend_length);

  printf("#file1: %s\n#file2: %s\n", peakset[0].file, peakset[1].file);
  printf("#num1: %d\tnum2: %d\tnum1_overlap: %d (%.1f%%)\tnum1_notoverlap: %d (%.1f%%)\tnum2_overlap: %d (%.1f%%)\tnum2_notoverlap: %d (%.1f%%)\n", peakset[0].num, peakset[1].num, peakset[0].cnt_overlap[1], (100*peakset[0].cnt_overlap[1]/(double)peakset[0].num), peakset[0].cnt_notoverlap[1], (100*peakset[0].cnt_notoverlap[1]/(double)peakset[0].num), peakset[1].cnt_overlap[0], (100*peakset[1].cnt_overlap[0]/(double)peakset[1].num), peakset[1].cnt_notoverlap[0], (100*peakset[1].cnt_notoverlap[0]/(double)peakset[1].num));
  printf("#peakwidth total1: %d bp\tpeakwidth total2: %d bp\toverlappeaks total: %d bp (%.2f%% / %.2f%%)\n", peakset[0].peakwid_total, peakset[1].peakwid_total, peakset[0].base_overlap[1], (peakset[0].base_overlap[1]/(double)peakset[0].peakwid_total)*100, (peakset[0].base_overlap[1]/(double)peakset[1].peakwid_total)*100);
  if(nobs) exit(0);

  if(unison==2){
    if(!showallcols) printf("chromosome\tstart\tend\tsummit 1\tsummit 2\tenrich 1\tenrich 2\tintensity 1\tintensity 2\tA (log2(1*2)/2)\tM (log2(1/2))\tdiff of summit\n");
    for(i=0; i<peakset[0].cnt_overlap_red[1]; i++){
      print1bs_overlap(peakset[0].bsarray_overlap[1], peakset[0].bsarray, peakset[1].bsarray, i);
    }
  }else if(unison==1){
    if(!showallcols) printf("chromosome\tstart\tend\tsummit\tlength\tenrich\tintensity\n");
    for(i=0; i<peakset[0].num; i++){
      if(peakset[0].bsarray[i].overlap[1]){
	print1bs(peakset[0].bsarray, i);
      }
    }
  }else{
    if(!showallcols) printf("chromosome\tstart\tend\tsummit\tlength\tenrich\tintensity\n");
    for(i=0; i<peakset[0].num; i++){
      if(!(peakset[0].bsarray[i].overlap[1])) print1bs(peakset[0].bsarray, i);
    }
  }

  for(i=0;i<samplenum;i++){
    free(peakset[i].bsarray);
    for(j=0; j<samplenum; j++) free(peakset[i].bsarray_overlap[j]);
    free(peakset[i].bsarray_overlap);
    free(peakset[i].base_overlap);
    free(peakset[i].cnt_overlap);
    free(peakset[i].cnt_notoverlap);
    free(peakset[i].cnt_overlap_red);
  }
  return 0;
}

static void argv_init(int argc, char **argv, struct peakset *peakset, int *extend_length, int samplenum){
  int i;
  unison=0;
  maxposi=0;
  nobs=0;
  sp_scer=0;
  showallcols=0;
  for(i=0; i<samplenum; i++){
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

  for(i=1; i<argc; i++){
    if(!strcmp(argv[i], "-1"))              peakset[0].file = argv[++i];
    else if(!strcmp(argv[i], "-2"))         peakset[1].file = argv[++i];
    else if(!strcmp(argv[i], "-n1"))        peakset[0].num_defined = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-n2"))        peakset[1].num_defined = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-l"))         *extend_length=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-red"))       unison=2;
    else if(!strcmp(argv[i], "-and"))       unison=1;
    else if(!strcmp(argv[i], "-not"))       unison=0;
    else if(!strcmp(argv[i], "-maxposi"))   maxposi=1;
    else if(!strcmp(argv[i], "-scer"))      sp_scer=1;
    else if(!strcmp(argv[i], "-nobs"))      nobs=1;
    else if(!strcmp(argv[i], "-showallcols"))      showallcols=1;
    else goto err;
  }
  if(!(peakset[0].file) || !(peakset[1].file)) goto err;
  return;

 err:
  fprintf(stderr, "%s", Usage);
  exit(0);
}
