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
Usage    :
    ./CLASSIFY_ONLY.sh [OPTION]

Trio phase TGS reads use exist kmer datasets.
Options  :
        --paternal_mer paternal unique kmers
        --maternal_mer maternal unique kmers
        --offspring    offspring TGS reads file in FASTA format.
                       file in gzip format can be accepted, but filename must end by .gz.
        --format       fasta/fastq . set the format of --offspring.
                       [ optional, default fasta. ]
        --thread       thread num.
                       [ optional, default 8 threads. ]
        --help         print this usage message.

Examples :
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --offspring son.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --offspring son.fastq --format fastq
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --offspring son.fasta --offspring son2.fasta

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --offspring son.fasta --thread 20
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
        --offspring   offspring TGS reads file in FASTA format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --format      fasta/fastq . set the format of --offspring.
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
    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq --offspring son.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq --offspring son.L01.fasta --offspring son.L02.fasta

    ./HAST4TGS.sh --paternal father.fastq --maternal mother.fastq \
                     --offspring son.fasta --memory 50 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=32 --m-lower=8 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux

```
### format of phased.out 

Below is the detail format of phased.out:
```
Read    read_full_name  classify_result read_length kmer_probaility_1   kmer_probaility_2  kmer_density_1  kmer_density_2  kmer_count_1 kmer_count_1
```

**If you want to parse this file by yourself, please double check whether your read_full_name contains multi-columns or not.**

```
# in HAST4TGS , I use paternal-specific kmer-lirary as lib1
kmer_probaility_1 = kmer_count_1 / kmer-number-in-lib1 
# in HAST4TGS , I use maternal-specific kmer-lirary as lib1
kmer_probaility_2 = kmer_count_2 / kmer-number-in-lib2 

kmer_density_1    = kmer_count_1 / read_length
kmer_density_2    = kmer_count_2 / read_length
```

The Pseudo code of detecting classify_result
```
if kmer_probaility_1 > kmer_probaility_2 
    classify_result = haplotype0
elseif kmer_probaility_1 < kmer_probaility_2
    classify_result = haplotype1
else
    classify_result = ambiguous
endif
```

Here is an example of phased.out:
```
Read    m64064_201214_094226/9/ccs      ambiguous       15723   0       0       0       0       0       0
Read    m64064_201214_094226/10/ccs     haplotype1      14475   2.01279e-07     2.4242e-05      0.00332065      0.233829        483380
Read    m64064_201214_094226/14/ccs     haplotype0      13666   2.60321e-05     0       0.454932        0       6208    0
Read    m64064_201214_094226/15/ccs     haplotype1      13572   0       3.22749e-07     0       0.00332054      0       45
Read    m64064_201214_094226/24/ccs     haplotype1      14796   0       1.17624e-06     0       0.0110991       0       164
Read    m64064_201214_094226/25/ccs     haplotype0      13647   2.51599e-08     7.17219e-09     0.000440302     7.33837e-05     6 1
Read    m64064_201214_094226/30/ccs     haplotype0      15676   2.9261e-05      0       0.445708        0       6978    0
```

Enjoy !
