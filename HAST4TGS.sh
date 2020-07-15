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
        --paternal    paternal NGS reads file in FASTQ format.
                      ( note : gzip format is NOT supported. )
        --maternal    maternal NGS reads file in FASTQ format.
                      ( note : gzip format is NOT supported. )
        --filial      filial TGS reads file in FASTA format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --format      fasta/fastq . set the format of --filial.
                      [ optional, default fasta. ]
        --thread      thread num.
                      [ optional, default 8 threads. ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      ( note: real memory used may be greater than this. )
                      [ optional, default 20GB. ]
        --jellyfish   jellyfish path.
                      [ optional, default jellyfish. ]
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
    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq --format fastq

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.L01.fasta --filial son.L02.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq \
                     --filial son.fasta --memory 20 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=32 --m-lower=8 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux
"""
}

###############################################################################
# basic variables 
###############################################################################
MER=21
JELLY=jellyfish
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
        "--jellyfish")
            JELLY=$2
            shift
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
            PATERNAL=$2
            shift
            ;;
        "--maternal")
            MATERNAL=$2
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
echo "    jellyfish      : $JELLY"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower(maternal): $MLOWER"
echo "    upper(maternal): $MUPPER"
echo "    lower(paternal): $PLOWER"
echo "    upper(paternal): $PUPPER"
echo "    auto_bounds    : $AUTO_BOUNDS"
echo "HAST.sh in dir  : $SPATH"

CLASSIFY=$SPATH"/classify"
ANALYSIS=$SPATH"/analysis_kmercount.sh"

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL || -z $FILIAL || \
    -z $JELLY  || $MER -lt 11 || \
    $MLOWER -lt 1 || $MUPPER -gt 100000000 || \
    $PLOWER -lt 1 || $PUPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
if [[ $FORMAT != 'fasta' && $FORMAT != 'fastq' ]] ; then 
    echo "ERROR : format invalid ... exit!!!"
    exit 1
fi
if [[ ! -e $CLASSIFY ]] ; then 
    echo "ERROR : please run \"make\" command in $SPATH before using this script! exit..."
    exit 1
fi
if [[ ! -e $ANALYSIS ]] ; then
    echo "ERROR : \"$ANALYSIS\"  is missing. please download it from github. exit..."
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
# count NGS reads
echo "extract unique mers by jellyfish ..."
if [[ ! -e "step_01_done" ]] ; then
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  maternal_mer_counts.jf $MATERNAL || exit 1
    date >> "step_01_done"
fi
if [[ ! -e "step_02_done" ]] ; then
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  paternal_mer_counts.jf $PATERNAL || exit 1
    date >> "step_02_done"
fi
# dump all mers
if [[ ! -e "step_03_done" ]] ; then
    $JELLY dump maternal_mer_counts.jf            -o maternal.mer.fa || exit 1
    date >> "step_03_done"
fi

if [[ ! -e "step_04_done" ]] ; then
    $JELLY dump paternal_mer_counts.jf            -o paternal.mer.fa || exit 1
    date >> "step_04_done"
fi

if [[ $AUTO_BOUNDS == 1 ]] ; then 
    if [[ ! -e "step_04_1_done" ]] ; then
        sh $ANALYSIS || exit 1
        date >> "step_04_1_done"
    fi
    MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
    PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
fi

echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "
# dump filter mers
if [[ ! -e "step_05_done" ]] ; then
    $JELLY dump -L $MLOWER -U $MUPPER maternal_mer_counts.jf -o maternal.mer.filter.fa || exit 1
    date >> "step_05_done"
fi
if [[ ! -e "step_06_done" ]] ; then
    $JELLY dump -L $PLOWER -U $PUPPER paternal_mer_counts.jf -o paternal.mer.filter.fa || exit 1
    date >> "step_06_done"
fi
if [[ ! -e "step_07_done" ]] ; then
    # rm temporary files
    rm maternal_mer_counts.jf paternal_mer_counts.jf
    date >> "step_07_done"
fi

if [[ ! -e "step_08_done" ]] ; then
    # rm temporary files
    # mix 1 copy of paternal mers and 2 copy of maternal mers
    cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa || exit 1
    rm maternal.mer.fa
    rm paternal.mer.fa
    date >> "step_08_done"
fi

if [[ ! -e 'step_09_done' ]] ; then
    # count p/maternal mixed mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o mixed_mer_counts.js mixed.fa || exit 1
    date >> 'step_09_done'
fi
if [[ ! -e "step_10_done" ]] ; then
    # count==1 refer to paternal unique mers
    $JELLY dump -U 1 mixed_mer_counts.js          >paternal.mer.unique.fa || exit 1
    date >> "step_10_done"
fi
if [[ ! -e "step_11_done" ]] ; then
    # count==2 refer to maternal unique mers
    $JELLY dump -L 2 -U 2 mixed_mer_counts.js     >maternal.mer.unique.fa || exit 1
    date >> "step_11_done"
fi
if [[ ! -e "step_12_done" ]] ; then
    # rm temporary files
    rm mixed.fa mixed_mer_counts.js
    date >> "step_12_done"
fi
if [[ ! -e "step_13_done" ]] ; then
    # mix unique mers and filter mers
    cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa || exit 1
    cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa || exit 1
    rm paternal.mer.unique.fa
    rm paternal.mer.filter.fa
    rm maternal.mer.filter.fa
    rm maternal.mer.unique.fa
    date >> "step_13_done"
fi
if [[ ! -e "step_14_done" ]] ; then
    # count unique and filer mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa || exit 1
    date >> "step_14_done"
fi
if [[ ! -e "step_15_done" ]] ; then
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa || exit 1
    date >> "step_15_done"
fi
if [[ ! -e "step_16_done" ]] ; then
    # extrat both unique and filter mers
    $JELLY dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.mer || exit 1
    date >> "step_16_done"
fi
if [[ ! -e "step_17_done" ]] ; then
    $JELLY dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.mer || exit 1
    date >> "step_17_done"
fi
if [[ ! -e "step_18_done" ]] ; then
    # rm temporary files
    rm paternal_mixed.mer.fa paternal_mixed_mer_counts.js
    rm maternal_mixed.mer.fa maternal_mixed_mer_counts.js
    date >> "step_18_done"
fi
echo "extract unique mers done..."
date
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
echo "extract unique barcode by classify ..."
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done
if [[ ! -e "step_19_done" ]] ; then
    $CLASSIFY --hap paternal.unique.filter.mer --hap maternal.unique.filter.mer \
        --thread $CPU  $READ --format $FORMAT >phasing.out  2>phasing.log || exit 1
    date >> "step_19_done"
fi

if [[ ! -e "step_20_done" ]] ; then
    cat phasing.out |grep Read |grep haplotype0 |awk '{print $2}' > paternal.cut
    cat phasing.out |grep Read |grep haplotype1 |awk '{print $2}' > maternal.cut
    cat phasing.out |grep Read |grep ambiguous |awk '{print $2}' > ambiguous.cut
    date >> "step_20_done"
fi
if [[ ! -e "step_21_done" ]] ; then
    # extract the reads
    if [[ $FORMAT == 'fasta' ]] ; then
        for x in $FILIAL
        do
            name=`basename $x`
            if [[ ${name: -3} == ".gz" ]] ; then
                name=${name%%.gz}
                gzip -dc $x | awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal.cut  - >"paternal."$name
                gzip -dc $x | awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal.cut  - >"maternal."$name
                gzip -dc $x | awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' ambiguous.cut - >"ambiguous."$name
            else 
                awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal.cut  $x >"paternal."$name
                awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal.cut  $x >"maternal."$name
                awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' ambiguous.cut $x >"ambiguous."$name
            fi
        done
    else
        for x in $FILIAL
        do
            name=`basename $x`
            if [[ ${name: -3} == ".gz" ]] ; then
                name=${name%%.gz}
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal.cut  - >"paternal."$name
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal.cut  - >"maternal."$name
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' ambiguous.cut - >"ambiguous."$name
            else 
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal.cut  $x >"paternal."$name
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal.cut  $x >"maternal."$name
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' ambiguous.cut $x >"ambiguous."$name
            fi
        done
    fi
    date >> "step_21_done"
fi
echo "phase reads done"
date
echo "__END__"
