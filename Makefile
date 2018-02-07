CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPCMNDIR = $(SSPDIR)/common
SSPOBJDIR = $(SSPDIR)/obj
SSPCMNOBJDIR = $(SSPDIR)/cobj
ALGLIBDIR = $(SSPCMNDIR)/alglib

PROGRAMS = compare_bed2loop gtf2refFlat compare_bed2tss peak_occurance multibed2gene compare_bs FRiR
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SSPCMNOBJDIR)/util.o $(SSPCMNOBJDIR)/ReadAnnotation.o
OBJS_GTF = $(OBJDIR)/gtf2refFlat.o
OBJS_LOOP = $(OBJDIR)/compare_bed2loop.o $(OBJDIR)/gene_bed.o
OBJS_COM = $(OBJDIR)/compare_bed2tss.o $(OBJDIR)/gene_bed.o
OBJS_PO = $(OBJDIR)/peak_occurance.o $(OBJDIR)/gene_bed.o $(SSPCMNOBJDIR)/statistics.o #$(ALGLIBDIR)/libalglib.a
OBJS_MG = $(OBJDIR)/multibed2gene.o $(OBJDIR)/gene_bed.o
OBJS_FRIR = $(OBJDIR)/FRiR.o $(OBJS_SSP) $(SSPOBJDIR)/Mapfile.o $(SSPOBJDIR)/LibraryComplexity.o $(SSPOBJDIR)/ParseMapfile.o $(SSPCMNOBJDIR)/BoostOptions.o $(SSPCMNOBJDIR)/util.o $(SSPOBJDIR)/ReadBpStatus.o 

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/gtf2refFlat: $(OBJS_GTF) $(OBJS_UTIL)
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
$(BINDIR)/compare_bs: src/compare_bs.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	gcc -o $@ $^ -Wall -W -O3 -lm
$(BINDIR)/FRiR: $(OBJS_FRIR)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS)


$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) $(OBJDIR)/$(notdir $@)
$(SSPCMNOBJDIR)/%.o: $(SSPCMNDIR)/%.cpp
	$(MAKE) -C $(SSPDIR) cobj/$(notdir $@)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf $(BINDIR) $(OBJDIR)

HEADS_UTIL = $(SSPCMNDIR)/util.hpp $(SSPCMNDIR)/inline.hpp $(SSPCMNDIR)/seq.hpp $(SSPCMNDIR)/BedFormat.hpp

$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_GTF): Makefile $(HEADS_UTIL)
$(OBJS_COM): Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_LOOP): Makefile $(HEADS_UTIL)
$(OBJS_PO):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_MG):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
