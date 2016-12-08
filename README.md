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
    git clone https://github.com/rnakato/ChIPseqTools.git
    git clone https://github.com/rnakato/SSP.git ChIPseqTools/src/SSP
    cd ChIPseqTools
    make

#### 2.4. Add the PATH environment variable
For example, if you downloaded DROMPA and cpdf into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/ChIPseqTools/bin

#3. Usage

#### 3.1. gtf2refFlat

#### 3.2. compare_bed2tss

#### 3.3. peak_occurance

#### 3.4. multibed2gene
