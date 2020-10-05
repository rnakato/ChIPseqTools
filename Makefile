.PHONY: all clean

BINDIR = bin
PROGRAMS_C = compare_bs mergebed2CRM
PROGRAMS_CPP = compare_bed2loop gtf2refFlat parseChIADropReadList compare_bed2tss peak_occurance multibed2gene FRiR

TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS_C) $(BINDIR)/,$(PROGRAMS_CPP))

SSPDIR = DROMPAplus/submodules/SSP
HTSLIBDIR = $(SSPDIR)/src/htslib-1.10.2/

all: $(TARGET) $(HTSLIBDIR)/libhts.a

$(TARGET): $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
#	cp build/test/parse2wig/parse2wig+ build/test/drompa/drompa+ bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
	make -C $(SSPDIR) clean
