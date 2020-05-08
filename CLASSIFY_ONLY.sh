#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """Options  :
        --paternal_mer paternal unique kmers
        --maternal_mer maternal unique kmers
        --filial       filial TGS reads file in FASTA format.
                       file in gzip format can be accepted, but filename must end by ".gz".
        --format       fasta/fastq . set the format of --filial.
                       [ optional, default fasta. ]
        --thread       thread num.
                       [ optional, default 8 threads. ]
        --help         print this usage message.

Examples :
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --format fastq
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta --filial son2.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta --thread 20
"""
}

###############################################################################
# basic variables 
###############################################################################
CPU=8
PATERNAL=""
MATERNAL=""
FILIAL=""
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
        "--thread")
            CPU=$2
            shift
            ;;
        "--paternal_mer")
            PATERNAL=$2
            shift
            ;;
        "--maternal_mer")
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
echo "CLASSIFY_ONLY starting with : "
echo "    paternal kmers : $PATERNAL"
echo "    maternal kmers : $MATERNAL"
echo "    filial input   : $FILIAL"
echo "    filial format  : $FORMAT"
echo "    thread         : $CPU "
echo "CLASSIFY_ONLY.sh in dir  : $SPATH"

CLASSIFY=$SPATH"/classify"

# sanity check
if [[ $CPU -lt 1 || -z $PATERNAL || -z $MATERNAL || -z $FILIAL ]] ; then
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
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
echo "extract unique barcode by classify ..."
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done
if [[ ! -e "step_19_done" ]] ; then
    $CLASSIFY --hap $PATERNAL --hap $MATERNAL \
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
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' paternal.cut  - >"paternal."$name
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' maternal.cut  - >"maternal."$name
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' ambiguous.cut - >"ambiguous."$name
            else 
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' paternal.cut  $x >"paternal."$name
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' maternal.cut  $x >"maternal."$name
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' ambiguous.cut $x >"ambiguous."$name
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
