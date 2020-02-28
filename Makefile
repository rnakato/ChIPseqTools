CC = g++
CFLAGS = -std=c++11 -O2 -Wall -W
CFLAGS_C  = -Wall -g -W -O3 -lm
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas

SRCDIR = ./src
SRC_CDIR = ./src_C
OBJDIR = ./obj
OBJ_CDIR = ./obj_C
LIBDIR = ./lib
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj
ALGLIBDIR = $(SSPCMNDIR)/alglib

PROGRAMS = compare_bed2loop gtf2refFlat parseChIADropReadList compare_bed2tss peak_occurance multibed2gene FRiR compare_bs mergebed2CRM
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/ReadAnnotation.o
OBJS_GTF  = $(OBJDIR)/gtf2refFlat.o
OBJS_PCD  = $(OBJDIR)/parseChIADropReadList.o
OBJS_LOOP = $(OBJDIR)/compare_bed2loop.o $(OBJDIR)/gene_bed.o
OBJS_COM  = $(OBJDIR)/compare_bed2tss.o $(OBJDIR)/gene_bed.o
OBJS_PO   = $(OBJDIR)/peak_occurance.o $(OBJDIR)/gene_bed.o $(SSPCMNOBJDIR)/statistics.o
OBJS_MG   = $(OBJDIR)/multibed2gene.o $(OBJDIR)/gene_bed.o
OBJS_FRIR = $(OBJDIR)/FRiR.o $(OBJS_SSP) $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/LibraryComplexity.o $(SSPOBJDIR)/ParseMapfile.o $(SSPCMNOBJDIR)/BoostOptions.o $(SSPCMNOBJDIR)/util.o $(SSPOBJDIR)/ReadBpStatus.o
OBJS_C_UTIL = $(OBJ_CDIR)/compare.o $(OBJ_CDIR)/readgene.o $(OBJ_CDIR)/my.o $(OBJ_CDIR)/stringp.o
OBJS_CB  = $(OBJ_CDIR)/compare_bs.o $(OBJS_C_UTIL)
OBJS_CRM = $(OBJ_CDIR)/mergebed2CRM.o $(OBJS_C_UTIL)

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/gtf2refFlat: $(OBJS_GTF) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/parseChIADropReadList: $(OBJS_PCD) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/compare_bed2loop: $(OBJS_LOOP) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/compare_bed2tss: $(OBJS_COM) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/peak_occurance: $(OBJS_PO) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/multibed2gene: $(OBJS_MG) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)
$(BINDIR)/FRiR: $(OBJS_FRIR)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)

### C programs
$(BINDIR)/compare_bs: $(OBJS_CB)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	gcc -o $@ $^ $(CFLAGS_C)
$(BINDIR)/mergebed2CRM: $(OBJS_CRM)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	gcc -o $@ $^ $(CFLAGS_C)

$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) $(OBJDIR)/$(notdir $@)
$(SSPCMNOBJDIR)/%.o: $(SSPCMNDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) cobj/$(notdir $@)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)
$(OBJ_CDIR)/%.o: $(SRC_CDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	gcc -o $@ -c $< $(CFLAGS_C)

clean:
	rm -rf $(BINDIR) $(OBJDIR) $(OBJ_CDIR)

HEADS_UTIL = $(SSPCMNDIR)/util.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/BedFormat.hpp

$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_GTF): Makefile $(HEADS_UTIL)
$(OBJS_COM): Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_LOOP): Makefile $(HEADS_UTIL)
$(OBJS_PO):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_MG):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
