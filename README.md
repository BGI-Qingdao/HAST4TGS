# HAST
A fast and memory effiective version of trio-binning .

- Use jellyfish to replace meryl 
- Rewrite classify.py by c++ and add multi-thread support .

## INSTALL

```
git clone https://github.com/BGI-Qingdao/HAST4TGS.git
cd HAST
make
```

## USAGE

```
Usage    :
Usage    :
    ./HAST4TGS.sh [OPTION]

A fast and memory effiective version of trio-binning .

Options  :
        --paternal    paternal NGS reads file in fastq format.
                      ( note : gzip format IS NOT supported. )
        --maternal    maternal NGS reads file in fastq format.
                      ( note : gzip format IS NOT supported. )
        --filial      filial stLFR reads file in fastq format.
                      file in gzip format is accepted, but filename must end by ".gz".
        --thread      threads num.
                      [ optional, default 8 thread. ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      ( note: real memory used maybe greater than this. )
                      [ optional, default 20GB. ]
        --jellyfish   jellyfish path.
                      [ optional, default jellyfish. ]
        --mer         mer-size
                      [ optional, default 21. ]
        --m-lower     maternal kmer count tablle will ignore mer with count < m-lower.
                      [ optional, default 9. ]
        --m-upper     maternal kmer count tablle will ignore mer with count > m-upper.
                      [ optional, default 33. ]
        --p-lower     paternal kmer count tablle will ignore mer with count < p-lower.
                      [ optional, default 9. ]
        --p-upper     paternal kmer count tablle will ignore mer with count > p-upper.
                      [ optional, default 33. ]
        --help        print this usage message.

Examples :
    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq --filial son.r1.fastq --filial son.r2.fastq

    ./HAST4TGS.sh --paternal father.fastq --maternal mater.fastq \
                     --filial son.r1.fastq --memory 20 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=33 --m-lower=9 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux

```

Enjoy !
