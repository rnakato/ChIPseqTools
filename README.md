# ChIPseqTools

#1. Overview
This repository contains several utility tools mainly for ChIP-seq analysis.
These programs are written in C++ using [Boost library](http://www.boost.org/).

#2. Install

#### 2.1. Install required libraries
for Ubuntu:

    sudo apt-get install git build-essential libboost-all-dev
 
for CentOS:

    sudo yum -y install git gcc-c++ boost-devel

#### 2.3. Install 
    git clone --recursive https://github.com/rnakato/ChIPseqTools.git
    cd ChIPseqTools
    make

#### 2.4. Add the PATH environment variable
For example, if you downloaded DROMPA and cpdf into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/ChIPseqTools/bin

#### 2.5 Update
 git pull origin master
 git submodule foreach git pull origin master

#3. Usage

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

#### 3.3. peak_occurance

#### 3.4. multibed2gene
