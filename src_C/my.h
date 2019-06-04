/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MY_H_
#define _MY_H_

#include <stdio.h>
#include <stdlib.h>

FILE *my_fopen_r(char *filename);
FILE *my_fopen_w(char *filename);
FILE *my_fopen_wb(char *filename);
FILE *my_fopen_a(char *filename);
void *my_calloc(size_t n, size_t s, char *name);
void *my_realloc(void *p, size_t s, char *name);
void my_system(char *);
void *inc_and_realloc(void *, int, int *, int, size_t, char *);

#define MYFREE(p) {if(p){free(p); (p)=NULL;} }

#endif
