/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "stringp.h"

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

int ParseLine_arbit(char *str, struct elem clm[], char token){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len, sizeof(char), "strtemp");
  for(i=0; i<=len; i++){
    if(str[i]=='\0'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      MYFREE(strtemp);
      return ++num;
    }
    if(str[i]==token){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++; 
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  MYFREE(strtemp);
  return num;
}

