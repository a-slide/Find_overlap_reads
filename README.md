# Find_overlap_reads

[see GitHub Page](http://a-slide.github.io/Find_overlap_reads) 

## Motivation
Find_overlap_reads is a **python3** object oriented script, developed to parse a genomic alignment files(BAM/SAM). The program is intended to count and extract reads or read pairs overlapping genomic intervals provided in a tab-separated text file.

## Principle

1. The interval file is parsed and Interval objects are created for each valid line within it.
2. Sam/Bam (CRAM?) files are parsed iteratively with *pysam*.
3. For each read the program verify if the read itself or the template covered by the read and its pair overlap one of the genomic interval.
4. A report containing the total number of reads, the number of read mapped and the reads overlapping each interval is created in the current dir.
5. A subset of bam containing only read overlapping intervals is generated, sorted and indexed.  

## Dependencies

* [pysam](https://github.com/pysam-developers/pysam) 0.8.1+ (based on htslib and samtools versions 1.1)

If you have pip already installed, enter the following line to install pysam:
```bash
sudo pip install pysam
```

## Get Find_overlap_reads

* Clone the repository with --recursive option to also pull the submodule
``` bash
$ git clone --recursive https://github.com/a-slide/Find_overlap_reads/ my_folder/
```

* Enter the root of the program folder and make the main script executable
``` bash
$ sudo chmod u+x Find_overlap_reads.py
```

* Add Find_overlap_reads.py in your PATH

## Usage

    Usage: find_overlap_read -f genomic_interval.csv [-b/-r] f1.bam(sam),[f2.bam(sam)...fn.bam(sam)
    Parse a BAM/SAM file(s) and extract reads overlapping given genomic coordinates

    Options:
        --version         show program's version number and exit
        -h, --help        show this help message and exit
        -f INTERVAL_FILE  Path of the tab separated file contaning genomic interval (mandatory)
        -t TEMPLATE_LEN   Maximum length of the template between 2 paired reads to be considered as a concordant pair (default = 1000)
        -b, --no_bam      Don't output bam file(s) (default = True)
        -r, --no_report   Don't output report file(s) (default = True)

#### Genomic interval file

The file containing genomic intervals have to be formated as a **tab-separated values** text format as follow:

```Name_of_the_the_reference    start_coordinate(INT)    end_coordinate(INT)    Name_of_the_interval (facultative)```

Lines with invalid integer values or less than 3 fields will be skipped. An error will be raised if no valid interval was found.
Examples of valid and invalid files are provided in the demo/ folder. A standard Bed file would be considered as a 

#### Aligment files (sam/bam)

A list of bam, sam and or cram files containing

## Development notebook

2 possibilities:
* Use ipython notebook with Dev_notebook.ipynb
* Consult directly online through nbviewer : [Notebook](http://nbviewer.ipython.org/github/a-slide/Find_overlap_reads/blob/master/Dev_notebook.ipynb)

## Authors and Contact

Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
