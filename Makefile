.PHONY: all clean

BINDIR = bin
PROGRAMS_C = compare_bs compare_bs2CRM
PROGRAMS_CPP = compare_bed2loop gtf2refFlat parseChIADropReadList compare_bed2tss peak_occurance multibed2gene FRiR

TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS_C) $(BINDIR)/,$(PROGRAMS_CPP))
TARGETFILE = $(addprefix build/src/C/,$(PROGRAMS_C)) $(addprefix build/src/C++/,$(PROGRAMS_CPP))

SSPDIR = submodules/DROMPAplus/submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2/

all: $(TARGET) $(HTSLIBDIR)/libhts.a

$(TARGET): $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp $(TARGETFILE) bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
	make -C $(SSPDIR) clean
