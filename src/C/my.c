/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "my.h"

FILE *my_fopen_r(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "r"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

FILE *my_fopen_w(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "w"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

FILE *my_fopen_wb(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "wb"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

FILE *my_fopen_a(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "a"))==NULL){
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

void my_system(char *command){
  int val = system(command);
  if(val) fprintf(stderr, "system: command error. %s\n", command);
}

void *inc_and_realloc(void *array, int num, int *nummax, int addnum, size_t size, char *str){
  if(num>=(*nummax)){
    (*nummax) += addnum;
    return my_realloc(array, (*nummax)*size, str);
  }
  else return array;
}
