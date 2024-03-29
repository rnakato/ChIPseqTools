# ChIPseqTools

## 1. Overview
This repository contains several utility tools for ChIP-seq and other epigenome analysis.
These programs are written in ANCI C and C++11 using [Boost library](http://www.boost.org/).

## 2. Install

### 2.1. Docker image

We recommend to use the latest Docker image from [DockerHub](https://hub.docker.com/r/rnakato/ssp_drompa).

#### 2.1.1. Docker
 To use docker command, type:

    docker pull rnakato/ssp_drompa
    docker run -it --rm rnakato/ssp_drompa <command>

**Note**: When using the docker image, it is necessary to mount the directory by ``-v`` option to access the input files as follows:

    docker run -it --rm -v $(pwd):/mnt rnakato/ssp_drompa parse2wig+ \
             -i /mnt/ChIP.bam -o ChIP --odir /mnt/parse2wigdir+ --gt /mnt/genometable.txt

This command mounts the current directory to ``/mnt`` directory in the container.
Please see also the [document of Docker](https://docs.docker.com/storage/volumes/).

#### 2.1.2. Singularity

Singularity can also be used to execute the docker image:

    singularity build ssp_drompa.sif docker://rnakato/ssp_drompa
    singularity exec ssp_drompa.sif <command>

Singularity mounts the current directory automatically. If you access the files in the other directory, please mount by `--bind` option, for instance:

    singularity exec --bind /work ssp_drompa.sif <command>

This command mounts `/work` directory.

### 2.2. Building from source

#### 2.2.1. Install required libraries
for Ubuntu:

    sudo apt install git build-essential libboost-all-dev libcurl4-gnutls-dev ibgtkmm-3.0-dev \
                     libgsl-dev liblzma-dev libz-dev libbz2-dev libgzstream0 libgzstream-dev cmake

for CentOS:

    sudo yum -y install git gcc-c++ boost-devel

On Mac:

     brew install gsl gtk gtkmm pkgconfig curl xz zlib boost cmake gzstream

#### 2.2.2. Install
    git clone --recursive https://github.com/rnakato/ChIPseqTools.git
    cd ChIPseqTools
    make

#### 2.2.3. Add the PATH environment variable

    export PATH = $PATH:(PATH_TO_ChIPseqTools)/ChIPseqTools/bin

#### 2.2.4 (Optional) Update repository

    git pull origin master
    git submodule foreach git pull origin master

# 3. Usage

#### 3.1. gtf2refFlat

#### 3.2. compare_bs

output shared peaks between bed1 and bed2

    compare_bs -1 <bed1> -2 <bed2> -and

output bed1-unique peaks

    compare_bs -1 <bed1> -2 <bed2> -not

output stats only

    compare_bs -1 <bed1> -2 <bed2> -and -nobs

include neighboring peaks within 5kbp as overlaped peaks

    compare_bs -1 <bed1> -2 <bed2> -and -l 5000

consider peak summit (default: whole peak region)

    compare_bs -1 <bed1> -2 <bed2> -and -maxposi

#### 3.2. compare_bed2tss

gene (gtf) and peaks

    compare_bed2tss -g <gtf> -b <peak> --gt <genome table>

gene (refFlat) and peaks

    compare_bed2tss -g <refFlat> --refFlat -b <peak> --gt <genome table>

comare with gene body

    compare_bed2tss -g <gtf> -b <peak> --gt <genome table> --mode 1

proportion of peaks against whole genome

    compare_bed2tss -g <gtf> -b <peak> --gt <genome table> --mode 2

proportion of peaks against whole genome (distinguish exon and instron)

    compare_bed2tss -g <gtf> -b <peak> --gt <genome table> --mode 2 --intron

#### 3.3. compare_bed2loop

    Usage: compare_bed2loop [option] --bed1 <1st bed> --bed2 <2nd bed> --loop <loop file> -o <output> -gt <genome_table>
    Options:
       --bed1 arg               1st bed file
       --bed2 arg               2nd bed file
       --loop arg               Loop file
       -l [ --length ] arg (=0) Extend length for overlap
       --gt arg                 Genome table (tab-delimited file describing the name
                                and length of each chromosome)
       --nobs                   do not output the overlapped loop list
       --hiccups                HICCups format as input (default: Mango)
       -h [ --help ]            print this message


#### 3.4. peak_occurance

#### 3.5. multibed2gene

#### 3.6. mergebed2CRM

      mergebed2CRM -i <bs file> -name <name> [-i <bs file> -name <name> ...]
         -l: extend length (default:0)
         -n: number of peaks for clustering (default:3000, setting 0 means use all peaks)
         -qnt: quantitative analysis

#### 3.6. FRiR

Repeat analysis.

       FRiR [option] -r <repeatfile> -i <inputfile> -o <output> --gt <genome_table>

<repeatfile> is the RepeatMasker file downloaded from the UCSC genome browser. FRiR can allow a gzipped repeat file.
