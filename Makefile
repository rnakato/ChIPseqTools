CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin
SSPDIR = ./src/SSP
SSPSRCDIR = $(SSPDIR)/src
SSPOBJDIR = $(SSPDIR)/obj
ALGLIBDIR = $(SSPSRCDIR)/alglib

PROGRAMS = gtf2refFlat compare_bed2tss peak_occurance multibed2gene compare_bs
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(SSPOBJDIR)/readdata.o $(SSPOBJDIR)/util.o
OBJS_GTF = $(OBJDIR)/gtf2refFlat.o
OBJS_COM = $(OBJDIR)/compare_bed2tss.o $(OBJDIR)/gene_bed.o
OBJS_PO = $(OBJDIR)/peak_occurance.o $(OBJDIR)/gene_bed.o $(ALGLIBDIR)/libalglib.a
OBJS_MG = $(OBJDIR)/multibed2gene.o $(OBJDIR)/gene_bed.o

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/gtf2refFlat: $(OBJS_GTF) $(OBJS_UTIL)
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

$(LIBDIR)/libalglib.a:
	make -C $(ALGLIBDIR)
$(SSPOBJDIR)/%.o: $(SSPSRCDIR)/%.cpp
	make -C $(SSPDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin obj

HEADS_UTIL = $(SSPSRCDIR)/util.h $(SSPSRCDIR)/readdata.h $(SSPSRCDIR)/macro.h $(SSPSRCDIR)/seq.h

$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_GTF): Makefile $(HEADS_UTIL)
$(OBJS_COM): Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_PO):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
$(OBJS_MG):  Makefile $(HEADS_UTIL) $(SRCDIR)/gene_bed.h
