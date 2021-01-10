#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./HAST4TGS.sh [OPTION]

A version of HAST for TGS long reads, fast and memory-effiective trio-binning.

Options  :
        --paternal    paternal NGS reads file in FASTA/Q format.
        --maternal    maternal NGS reads file in FASTA/Q format.
        --filial      filial TGS reads file in FASTA format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --format      fasta/fastq . set the format of --filial.
                      [ optional, default fasta. ]
        --thread      thread num.
                      [ optional, default 8 threads. ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      ( note: real memory used may be greater than this. )
                      [ optional, default 20GB. ]
        --mer         mer-size
                      [ optional, default 21. ]
        --m-lower     maternal kmer frequency table will ignore kmers with count < m-lower.
                      [ optional, default 9. ]
        --m-upper     maternal kmer frequency table will ignore kmers with count > m-upper.
                      [ optional, default 33. ]
        --p-lower     paternal kmer frequency table will ignore kmers with count < p-lower.
                      [ optional, default 9. ]
        --p-upper     paternal kmer frequency table will ignore kmers with count > p-upper.
                      [ optional, default 33. ]
        --auto_bounds automatically calcuate lower and upper bounds based on kmer analysis.
                      [ optional, default not trigger; no parameter. ]
                      ( note : if auto_bounds is on, it will overwrite --*-lower and --*-upper  ]
        --help        print this usage message.

Examples :
    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fasta --auto_bounds

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq --format fastq --auto_bounds

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.L01.fasta --filial son.L02.fasta --auto_bounds

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq \
                     --filial son.fasta --memory 20 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=32 --m-lower=8 --p-upper=33
"""
}

###############################################################################
# basic variables 
###############################################################################
MER=21
CPU=8
MEMORY=10
PLOWER=9
PUPPER=33
MLOWER=9
MUPPER=33
PATERNAL=""
MATERNAL=""
FILIAL=""
AUTO_BOUNDS=0
FORMAT='fasta'
SPATH=`dirname $0`
###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 ]] ; then 
    usage
    exit 0
fi
echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "-h")
            usage
            exit 0
            ;;
        "--help")
            usage
            exit 0
            ;;
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--m-lower")
            MLOWER=$2
            shift
            ;;
        "--m-upper")
            MUPPER=$2
            shift
            ;;
        "--p-lower")
            PLOWER=$2
            shift
            ;;
        "--p-upper")
            PUPPER=$2
            shift
            ;;
        "--mer")
            MER=$2
            shift
            ;;
        "--auto_bounds")
            AUTO_BOUNDS=1
            ;;
        "--paternal")
            PATERNAL=$2" "$PATERNAL
            shift
            ;;
        "--maternal")
            MATERNAL=$2" "$MATERNAL
            shift
            ;;
        "--filial")
            FILIAL=$2" "$FILIAL
            shift 
            ;;
        "--format")
            FORMAT=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done
# print arguments
echo "HAST starting with : "
echo "    paternal input : $PATERNAL"
echo "    maternal input : $MATERNAL"
echo "    filial input   : $FILIAL"
echo "    filial format  : $FORMAT"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower(maternal): $MLOWER"
echo "    upper(maternal): $MUPPER"
echo "    lower(paternal): $PLOWER"
echo "    upper(paternal): $PUPPER"
echo "    auto_bounds    : $AUTO_BOUNDS"
echo "HAST.sh in dir  : $SPATH"

STEP0=$SPATH"/00.build_unshare_kmers_by_meryl/build_unshared_kmers.sh"
STEP1=$SPATH"/01.classify_reads/CLASSIFY_ONLY.sh"
CLASSIFY=$SPATH"/01.classify_reads/classify"

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL || -z $FILIAL || \
    $MER -lt 11 || \
    $MLOWER -lt 1 || $MUPPER -gt 100000000 || \
    $PLOWER -lt 1 || $PUPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
if [[ $FORMAT != 'fasta' && $FORMAT != 'fastq' ]] ; then 
    echo "ERROR : format invalid ... exit!!!"
    exit 1
fi
if [[ ! -e $STEP0 || ! -e $STEP1 ]] ; then
    echo "downloading failed ? please re-download this!!!"
    exit 1
fi
if [[ ! -e $CLASSIFY ]] ; then 
    echo "ERROR : please run \"make\" command in $SPATH/01.classify_reads before using this script! exit..."
    exit 1
fi

for x in $MATERNAL $PATERNAL $FILIAL
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done
date
echo "__START__"


###############################################################################
# extract paternal.unique.filter.mer & maternal.unique.filter.mer
###############################################################################
echo "build unshared-kmer by meryl  ..."
date
if [[ $AUTO_BOUNDS==1 ]] ; then
    $STEP0 --paternal "$PATERNAL" \
           --maternal "$MATERNAL" \
           --thread   $CPU \
           --memory   $MEMORY \
           --auto_bounds >step0.log 2>step0.err
else
    $STEP0 --paternal "$PATERNAL" \
           --maternal "$MATERNAL" \
           --thread   $CPU    \
           --memory   $MEMORY \
           --m-lower  $MLOWER \
           --m-upper  $MUPPER \
           --p-lower  $PLOWER \
           --p-upper  $PUPPER >step0.log 2>step0.err
fi
date
echo "build unshared-kmer done"
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
echo "phase reads by classify ..."
date
$STEP1 --paternal_mer paternal.unique.filter.mer \
       --maternal_mer maternal.unique.filter.mer \
       --filial       $FILIAL \
       --format       $FORMAT \
       --thread       $CPU  >step1.log 2>step1.err
echo "phase reads done"
date
echo "__END__"
