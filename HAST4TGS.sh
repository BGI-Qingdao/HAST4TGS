#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./HAST4TGS.sh [OPTION]" 
    echo ""
    echo "A fast and memory effiective version of trio-binning ."
    echo ""
    echo "Options  :"
    echo "        --paternal    paternal NGS reads file in fastq format."
    echo "                      ( note : gzip format IS NOT supported. ) "
    echo "        --maternal    maternal NGS reads file in fastq format."
    echo "                      ( note : gzip format IS NOT supported. ) "
    echo "        --filial      filial stLFR reads file in fastq format."
    echo "                      file in gzip format is accepted, but filename must end by \".gz\"."
    echo "        --thread      threads num."
    echo "                      [ optional, default 8 thread. ]"
    echo "        --memory      x (GB) of memory to initial hash table by jellyfish."
    echo "                      ( note: real memory used maybe greater than this. )"
    echo "                      [ optional, default 20GB. ]"
    echo "        --jellyfish   jellyfish path."
    echo "                      [ optional, default jellyfish. ]"
    echo "        --mer         mer-size"
    echo "                      [ optional, default 21. ]"
    echo "        --m-lower     maternal kmer count tablle will ignore mer with count < m-lower."
    echo "                      [ optional, default 9. ]"
    echo "        --m-upper     maternal kmer count tablle will ignore mer with count > m-upper."
    echo "                      [ optional, default 33. ]"
    echo "        --p-lower     paternal kmer count tablle will ignore mer with count < p-lower."
    echo "                      [ optional, default 9. ]"
    echo "        --p-upper     paternal kmer count tablle will ignore mer with count > p-upper."
    echo "                      [ optional, default 33. ]"
    echo "        --help        print this usage message."
    echo "        "
    echo "Examples :"
    echo "    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq"
    echo ""
    echo "    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.r1.fastq --filial son.r2.fastq"
    echo ""
    echo "    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq \\"
    echo "                     --filial son.r1.fastq --memory 20 --thread 20 \\"
    echo "                     --mer 21 --p-lower=9 --p-upper=33 --m-lower=9 --p-upper=33 \\"
    echo "                     --jellyfish /home/software/jellyfish/jellyfish-linux"
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
echo "    jellyfish      : $JELLY"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower(maternal): $MLOWER"
echo "    upper(maternal): $MUPPER"
echo "    lower(paternal): $PLOWER"
echo "    upper(paternal): $PUPPER"
echo "HAST.sh in dir  : $SPATH"

CLASSIFY=$SPATH"/classify"

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL || -z $FILIAL || \
    -z $JELLY  || $MER -lt 11 || \
    $MLOWER -lt 1 || $MUPPER -gt 100000000 || \
    $PLOWER -lt 1 || $PUPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
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
# extract paternal.unique.filter.mer & maternal.unique.filter.mer
###############################################################################
# count NGS reads
echo "extract unique mers by jellyfish ..."
$JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  maternal_mer_counts.jf $MATERNAL
$JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  paternal_mer_counts.jf $PATERNAL
# dump all mers
$JELLY dump maternal_mer_counts.jf            -o maternal.mer.fa
$JELLY dump paternal_mer_counts.jf            -o paternal.mer.fa

# dump filter mers
$JELLY dump -L $MLOWER -U $MUPPER maternal_mer_counts.jf -o maternal.mer.filter.fa
$JELLY dump -L $PLOWER -U $PUPPER paternal_mer_counts.jf -o paternal.mer.filter.fa
# rm temporary files
rm maternal_mer_counts.jf paternal_mer_counts.jf
# mix 1 copy of paternal mers and 2 copy of maternal mers
cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa
# count p/maternal mixed mers
$JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o mixed_mer_counts.js mixed.fa
# count==1 refer to paternal unique mers
$JELLY dump -U 1 mixed_mer_counts.js          >paternal.mer.unique.fa
# count==2 refer to maternal unique mers
$JELLY dump -L 2 -U 2 mixed_mer_counts.js     >maternal.mer.unique.fa
# rm temporary files
rm mixed.fa mixed_mer_counts.js
# mix unique mers and filter mers
cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa
cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa
# count unique and filer mers
$JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa
$JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa
# extrat both unique and filter mers
$JELLY dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.mer
$JELLY dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.mer
# rm temporary files
rm paternal_mixed.mer.fa paternal_mixed_mer_counts.js
rm maternal_mixed.mer.fa maternal_mixed_mer_counts.js
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
$CLASSIFY --hap paternal.unique.filter.mer --hap maternal.unique.filter.mer \
    --thread $CPU  $READ >phasing.out  2>phasing.log

cat phasing.out |grep Read |grep haplotype0 |awk '{print $2}' > paternal.cut
cat phasing.out |grep Read |grep haplotype1 |awk '{print $2}' > maternal.cut
cat phasing.out |grep Read |grep ambiguous |awk '{print $2}' > ambiguous.cut
# extract the reads

for x in $FILIAL
do
    name=`basename $x`
    if [[ ${name: -3} == ".gz" ]] ; then 
        gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' paternal.cut  - >"maternal."$name
        gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' maternal.cut  - >"paternal."$name
        gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' ambiguous.cut - >"ambiguous."$name
    else 
        awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' paternal.cut  $x >"maternal."$name
        awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' maternal.cut  $x >"paternal."$name
        awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %2==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; c=0}  } } }' ambiguous.cut $x >"ambiguous."$name
    fi
done
echo "phase reads done"
date
echo "__END__"
