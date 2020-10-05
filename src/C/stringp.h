/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _STRINGP_H_
#define _STRINGP_H_

#include <string.h>
#include "my.h"

#define STR_LEN 100000
#define ELEM_NUM 256
#define ELEM_NUM2 2560

struct elem{
  char str[8000];
};

void chomp(char *);
int ParseLine(char *, struct elem clm[]);
int ParseLine_arbit(char *, struct elem clm[], char);

#endif
