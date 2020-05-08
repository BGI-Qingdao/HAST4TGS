# HAST4TGS
A version of HAST for TGS long reads, fast and memory-effiective trio-binning.

- Use jellyfish to generate and calculate k-mers instead of meryl; 
- Rewrite classify.py by c++ to speed up and add multi-thread support;
- Add a function to automatically generate k-mer frequency threholds. (default off) 

## INSTALL

```
git clone https://github.com/BGI-Qingdao/HAST4TGS.git
cd HAST4TGS
make
```

## USAGE

### Working with exist kmers datasets

```
Options  :
        --paternal_mer paternal unique kmers
        --maternal_mer maternal unique kmers
        --filial       filial TGS reads file in FASTA format.
                       file in gzip format can be accepted, but filename must end by .gz.
        --format       fasta/fastq . set the format of --filial.
                       [ optional, default fasta. ]
        --thread       thread num.
                       [ optional, default 8 threads. ]
        --memory       x (GB) of memory to initial hash table by jellyfish.
                       ( note: real memory used may be greater than this. )
                       [ optional, default 20GB. ]
        --help         print this usage message.

Examples :
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --format fastq
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta --filial son2.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fasta
                                        --memory 20 --thread 20
```

### Full pipeline

```
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
    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq --filial son.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq --filial son.L01.fasta --filial son.L02.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq \
                     --filial son.fasta --memory 50 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=32 --m-lower=8 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux

```

Enjoy !
