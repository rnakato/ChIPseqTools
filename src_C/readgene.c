#include "readgene.h"

static char *checkchrname(char *name)
{
  char *p = strstr(name, "chr");
  if (p) return p+3;
  else return name;
}

static void changechr_yeast(struct bs **bsarray, int num, char *str, int sp_scer)
{
  char *chrstr= checkchrname(str);
  if (!strcmp(chrstr, "I"))         strcpy((*bsarray)[num].chr, "1");
  else if (!strcmp(chrstr, "II"))   strcpy((*bsarray)[num].chr, "2");
  else if (!strcmp(chrstr, "III"))  strcpy((*bsarray)[num].chr, "3");
  else if (!strcmp(chrstr, "IV"))   strcpy((*bsarray)[num].chr, "4");
  else if (!strcmp(chrstr, "V"))    strcpy((*bsarray)[num].chr, "5");
  else if (!strcmp(chrstr, "VI"))   strcpy((*bsarray)[num].chr, "6");
  else if (!strcmp(chrstr, "VII"))  strcpy((*bsarray)[num].chr, "7");
  else if (!strcmp(chrstr, "VIII")) strcpy((*bsarray)[num].chr, "8");
  else if (!strcmp(chrstr, "IX"))   strcpy((*bsarray)[num].chr, "9");
  else if (!strcmp(chrstr, "X")){
    if (sp_scer) strcpy((*bsarray)[num].chr, "10"); else strcpy((*bsarray)[num].chr, "X");
  }else if (!strcmp(chrstr, "XI"))   strcpy((*bsarray)[num].chr, "11");
  else if (!strcmp(chrstr, "XII"))  strcpy((*bsarray)[num].chr, "12");
  else if (!strcmp(chrstr, "XIII")) strcpy((*bsarray)[num].chr, "13");
  else if (!strcmp(chrstr, "XIV"))  strcpy((*bsarray)[num].chr, "14");
  else if (!strcmp(chrstr, "XV"))   strcpy((*bsarray)[num].chr, "15");
  else if (!strcmp(chrstr, "XVI"))  strcpy((*bsarray)[num].chr, "16");
  else strcpy((*bsarray)[num].chr, chrstr);
  return;
}

int read_bs(char *filename, struct bs **bsarray, int *bsarraynum, int sp_scer)
{
  int num=0, linenum;
  char *chrstr=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  struct elem clm[ELEM_NUM];

  FILE *IN = my_fopen_r(filename);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if (str[0]=='\n' || str[0]=='#') continue;
    chomp(str);

    linenum = ParseLine(str, clm);
    if (linenum<3) continue;
    if (strstr(clm[0].str, "num1")) continue;
    if (!strcmp(clm[0].str, "chromosome")) continue;
    chrstr = checkchrname(clm[0].str);
    changechr_yeast(bsarray, num, chrstr, sp_scer);
    (*bsarray)[num].start = atoi(clm[1].str);
    (*bsarray)[num].end = atoi(clm[2].str);
    if (clm[3].str){
      (*bsarray)[num].maxposi = atoi(clm[3].str);
      if (atoi(clm[3].str)< (*bsarray)[num].start || atoi(clm[3].str)> (*bsarray)[num].end) (*bsarray)[num].maxposi = ((*bsarray)[num].start + (*bsarray)[num].end)/2;
    }
    if (linenum>=6) (*bsarray)[num].enrich  = atof(clm[5].str);
    if (linenum>=7) (*bsarray)[num].maxIP   = atof(clm[6].str);
    strcpy((*bsarray)[num].line, str);

    num++;
    if (num>=*bsarraynum){
      *bsarraynum += STRUCT_BS_MAX;
      *bsarray = (struct bs *)my_realloc(*bsarray, *bsarraynum * sizeof(struct bs), "bsarray");
    }
  }
  fclose(IN);
  free(str);
  return num;
}

int read_genelist(char *genefile, struct genename **genelist)
{
  FILE *IN = my_fopen_r(genefile);
  int num=0;
  char *str = (char *)my_calloc(STR_LEN, sizeof(struct elem), "str");
  struct elem clm[ELEM_NUM];

  while((fgets(str, STR_LEN, IN))!=NULL){
    if (str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    strcpy((*genelist)[num].name, clm[0].str);
    num++;
    if (num>=STRUCT_GENE_MAX){
      fprintf(stderr,"STRUCT_GENE_MAX over.\n"); exit(1);
    }
  }
  fclose(IN);
  free(str);
  return num;
}

int read_gene_ENS(char *filename, struct gene **ref_gene, int sp_scer)
{
  FILE *IN;
  int i, num=0;
  char *str, *chrstr=NULL;
  str = (char *)my_calloc(STR_LEN, sizeof(struct elem), "str");
  struct elem clm[ELEM_NUM];
  struct elem *clm2 = (struct elem *)my_calloc(ELEM_NUM2, sizeof(struct elem), "clm2");
  struct gene *gene = *ref_gene;
  gene = (struct gene *)my_calloc(STRUCT_GENE_MAX, sizeof(struct gene), "gene");

  IN = my_fopen_r(filename);
  while ((fgets(str, STR_LEN, IN))!=NULL) {
    if (str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);

    chrstr = checkchrname(clm[5].str);
    if (!strcmp(chrstr, "I"))         strcpy(gene[num].chr, "1");
    else if (!strcmp(chrstr, "II"))   strcpy(gene[num].chr, "2");
    else if (!strcmp(chrstr, "III"))  strcpy(gene[num].chr, "3");
    else if (!strcmp(chrstr, "IV"))   strcpy(gene[num].chr, "4");
    else if (!strcmp(chrstr, "V"))    strcpy(gene[num].chr, "5");
    else if (!strcmp(chrstr, "VI"))   strcpy(gene[num].chr, "6");
    else if (!strcmp(chrstr, "VII"))  strcpy(gene[num].chr, "7");
    else if (!strcmp(chrstr, "VIII")) strcpy(gene[num].chr, "8");
    else if (!strcmp(chrstr, "IX"))   strcpy(gene[num].chr, "9");
    else if (!strcmp(chrstr, "X")){
      if (sp_scer) strcpy(gene[num].chr, "10"); else strcpy(gene[num].chr, "X");
    }else if (!strcmp(chrstr, "XI"))   strcpy(gene[num].chr, "11");
    else if (!strcmp(chrstr, "XII"))  strcpy(gene[num].chr, "12");
    else if (!strcmp(chrstr, "XIII")) strcpy(gene[num].chr, "13");
    else if (!strcmp(chrstr, "XIV"))  strcpy(gene[num].chr, "14");
    else if (!strcmp(chrstr, "XV"))   strcpy(gene[num].chr, "15");
    else if (!strcmp(chrstr, "XVI"))  strcpy(gene[num].chr, "16");
    else strcpy(gene[num].chr, chrstr);

    strcpy(gene[num].name, clm[0].str);

    if (!strcmp(clm[1].str, "protein_coding"))  gene[num].genetype = CODINGGENE;
    else if (!strcmp(clm[1].str, "pseudogene")) gene[num].genetype = PSEUDO;
    else if (!strstr(clm[1].str, "RNA"))        gene[num].genetype = RNA;
    else gene[num].genetype = OTHERS;
    gene[num].dir = atoi(clm[2].str);
    gene[num].start = atoi(clm[3].str);
    gene[num].end   = atoi(clm[4].str);
    strcpy(gene[num].desc, clm[6].str);
    gene[num].exonnum = atoi(clm[7].str);
    gene[num].exon = (struct exon *)my_calloc(gene[num].exonnum, sizeof(struct exon), "exon");
    ParseLine_arbit(clm[8].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].start = atoi(clm2[i].str);
    ParseLine_arbit(clm[9].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].end = atoi(clm2[i].str);
    strcpy(gene[num].ID, clm[9].str);

    num++;
    if (num>=STRUCT_GENE_MAX) {
      fprintf(stderr,"STRUCT_GENE_MAX over.\n");
      exit(1);
    }
  }
  fclose(IN);
  free(str);
  free(clm2);

  *ref_gene = gene;
  return num;
}

int read_refFlat(char *filename, struct gene **ref_gene, int sp_scer)
{
  FILE *IN;
  int i, num=0, nummax=STRUCT_GENE_MAX;
  char *chrstr=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(struct elem), "str");
  struct elem clm[ELEM_NUM];
  struct elem *clm2 = (struct elem *)my_calloc(ELEM_NUM2, sizeof(struct elem), "clm2");
  struct gene *gene = *ref_gene;
  gene = (struct gene *)my_calloc(nummax, sizeof(struct gene), "gene");

  IN = my_fopen_r(filename);
  while ((fgets(str, STR_LEN, IN))!=NULL){
    if (str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);

    chrstr = checkchrname(clm[2].str);
    if (!strcmp(chrstr, "I"))         strcpy(gene[num].chr, "1");
    else if (!strcmp(chrstr, "II"))   strcpy(gene[num].chr, "2");
    else if (!strcmp(chrstr, "III"))  strcpy(gene[num].chr, "3");
    else if (!strcmp(chrstr, "IV"))   strcpy(gene[num].chr, "4");
    else if (!strcmp(chrstr, "V"))    strcpy(gene[num].chr, "5");
    else if (!strcmp(chrstr, "VI"))   strcpy(gene[num].chr, "6");
    else if (!strcmp(chrstr, "VII"))  strcpy(gene[num].chr, "7");
    else if (!strcmp(chrstr, "VIII")) strcpy(gene[num].chr, "8");
    else if (!strcmp(chrstr, "IX"))   strcpy(gene[num].chr, "9");
    else if (!strcmp(chrstr, "X")){
      if (sp_scer) strcpy(gene[num].chr, "10"); else strcpy(gene[num].chr, "X");
    }else if (!strcmp(chrstr, "XI"))   strcpy(gene[num].chr, "11");
    else if (!strcmp(chrstr, "XII"))  strcpy(gene[num].chr, "12");
    else if (!strcmp(chrstr, "XIII")) strcpy(gene[num].chr, "13");
    else if (!strcmp(chrstr, "XIV"))  strcpy(gene[num].chr, "14");
    else if (!strcmp(chrstr, "XV"))   strcpy(gene[num].chr, "15");
    else if (!strcmp(chrstr, "XVI"))  strcpy(gene[num].chr, "16");
    else strcpy(gene[num].chr, chrstr);

    if (strcmp(clm[0].str, "")){
      strcpy(gene[num].name, clm[0].str);
    }else strcpy(gene[num].name, clm[1].str);
    strcpy(gene[num].ID, clm[1].str);
    if (strstr(clm[1].str, "NM")) strcpy(gene[num].type, "mRNA"); else strcpy(gene[num].type, "non-coding RNA");
    if (!strcmp(clm[3].str,"+")) gene[num].dir = 1; else gene[num].dir = -1;
    gene[num].start = atoi(clm[4].str);
    gene[num].end   = atoi(clm[5].str);
    //    strcpy(gene[num].desc, clm[6].str);
    gene[num].exonnum = atoi(clm[8].str);
    gene[num].exon = (struct exon *)my_calloc(gene[num].exonnum, sizeof(struct exon), "exon");
    ParseLine_arbit(clm[9].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].start = atoi(clm2[i].str);
    ParseLine_arbit(clm[10].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].end = atoi(clm2[i].str);
    num++;
    if (num>=nummax) {
      nummax += STRUCT_GENE_MAX;
      gene = (struct gene *)my_realloc(gene, nummax*sizeof(struct gene), "gene");
    }
  }
  fclose(IN);
  free(str);
  free(clm2);
  *ref_gene = gene;
  return num;
}

int read_refseq(char *filename, struct gene **ref_gene)
{
  FILE *IN;
  int i, num=0;
  char *str = (char *)my_calloc(STR_LEN, sizeof(struct elem), "str");
  struct elem clm[ELEM_NUM];
  struct elem *clm2 = (struct elem *)my_calloc(ELEM_NUM2, sizeof(struct elem), "clm2");
  struct gene *gene = *ref_gene;
  gene = (struct gene *)my_calloc(STRUCT_GENE_MAX, sizeof(struct gene), "gene");

  IN = my_fopen_r(filename);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if (str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);

    strcpy(gene[num].chr, checkchrname(clm[5].str));
    if (strcmp(clm[0].str, "")) strcpy(gene[num].name, clm[0].str);
    else strcpy(gene[num].name, clm[10].str);
    strcpy(gene[num].type, clm[1].str);
    if (!strcmp(clm[2].str,"+")) gene[num].dir = 1;
    else gene[num].dir = -1;
    gene[num].start = atoi(clm[3].str);
    gene[num].end   = atoi(clm[4].str);
    strcpy(gene[num].desc, clm[6].str);
    gene[num].exonnum = atoi(clm[7].str);
    gene[num].exon = (struct exon *)my_calloc(gene[num].exonnum, sizeof(struct exon), "exon");
    ParseLine_arbit(clm[8].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].start = atoi(clm2[i].str);
    ParseLine_arbit(clm[9].str, clm2, ',');
    for (i=0; i<gene[num].exonnum; i++) gene[num].exon[i].end = atoi(clm2[i].str);
    num++;
    if (num>=STRUCT_GENE_MAX){fprintf(stderr,"STRUCT_GENE_MAX over.\n"); exit(1);}
  }
  fclose(IN);
  free(str);
  free(clm2);
  *ref_gene = gene;
  return num;
}

void print1bs(struct bs *bsarray, int i)
{
  printf("chr%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",
         bsarray[i].chr, bsarray[i].start, bsarray[i].end,
         bsarray[i].maxposi,
         bsarray[i].end-bsarray[i].start,
         bsarray[i].enrich,
         bsarray[i].maxIP);
}

void print1bsstr(struct bs *bsarray, int i)
{
  printf("%s\n", bsarray[i].line);
}


void print1bs_overlap(struct bs_overlap *bsarray, struct bs *bs1, struct bs *bs2, int i)
{
  int s1 = bsarray[i].bs1;
  int s2 = bsarray[i].bs2;
  printf("chr%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", bs1[s1].chr, bs1[s1].start, bs1[s1].end, bs1[s1].maxposi, bs2[s2].maxposi, bs1[s1].enrich, bs2[s2].enrich, bs1[s1].maxIP, bs2[s2].maxIP, log2(bs1[s1].maxIP*bs2[s2].maxIP)/2, log2(bsarray[i].ratio_maxIP), bsarray[i].peaksummitdiff);
}

void print1bsstr_overlap(struct bs_overlap *bsarray, struct bs *bs1, struct bs *bs2, int i)
{
  int s1 = bsarray[i].bs1;
  int s2 = bsarray[i].bs2;
  printf("%s\t%s\n", bs1[s1].line, bs2[s2].line);
}
