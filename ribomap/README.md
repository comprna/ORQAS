Overview
------
Ribomap is a package that generates isoform-level ribosome profiles from ribosome profiling data. Ribosome profiling is a recently developed high-throughput sequencing technique that captures approximately 30 bp long ribosome-protected mRNA fragments during translation. Because of alternative splicing and genomic repetitive sequences, a ribosome-protected read may map to many places in the transcriptome, leading to discarded or arbitrary mappings when standard approaches are used. Ribomap addresses this problem by assigning reads to potential origins in the transcriptome proportional to the estimated transcript abundance. This results in a more accurate estimation of the ribosome pileup compared to naive read assignments.

Prerequisites for Ribomap
------
<!---
* [__FASTX-Toolkit__] (http://hannonlab.cshl.edu/fastx_toolkit/index.html) for preprocessing reads
-->
* [__Sailfish__ (latest improved version: Salmon v0.2.7)](https://github.com/kingsfordgroup/sailfish/releases/tag/v0.2.7) for transcript abundance estimation
* [__STAR__ (v2.4.0j)](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0j) for read mapping

You can use `include_prerequisites.sh` in the `scripts` folder to download the pre-compiled Sailfish and STAR executables from the github page. 

    cd scripts
    ./include_prerequisites.sh os_type

Where `os_type` can be either `linux` or `osx`.

If you download the Ribomap binaries from the [release page](https://github.com/Kingsford-Group/ribomap/releases/tag/v1.0), STAR and Salmon executables and libraries are already included in the tar ball. 

Please make sure that the `PATH` variabile contains `<ribomap_dir>/bin`, and that `LD_LIBRARY_PATH` (or `DYLD_FALLBACK_LIBRARY_PATH` on OSX) contains `<ribomap_dir>/lib`.

NOTE: The pre-compiled executables might not work for all versions of operating systems. If Salmon takes an extremely long time to finish, please refer to the [online document](http://sailfish.readthedocs.org/en/develop/salmon.html) of Salmon and build Salmon from source. 

Compile from Source code
------
### Prerequisites
* [boost](http://www.boost.org/)
* [seqan (v1.4.1)](http://www.seqan.de/)

### Compile
a C++ compiler that support c++11 features (for instance g++ >= 4.7) is required.

    cd src
    make riboprof INC="-I/opt/local/include"
    make install

This will generate a c++ executable `riboprof` that assign ribosome profiling reads to transcript locations and copy the executable to the `bin` directory.

Please add the path for the prerequistite headers with flag `INC="-I<path/to/include/>"`

Run Ribomap
------
#### Run Ribomap with automatic transcript abundance estimation
`run_ribomap.sh` is an automatic pipeline for ribosome profiling data. It takes in the riboseq data and the RNA-seq data and automatically estimates the transcript abundance, then assigns riboseq reads to transcript locations based on the estimated transcript abundance. 

Under the `scripts` directory, run:

      ./run_ribomap.sh [options]

The list of options are as follows:
* __--rnaseq_fq__ (required) Input RNA-seq read fastq.gz file for transcript abundance estimation.
* __--riboseq_fq__ (required) Input ribosome profiling (riboseq) read fastq.gz file.
* __--transcript_fa__ (required) Input trascriptome reference fasta file.
* __--contaminant_fa__ Input contaminant sequence fasta file.
* __--cds_range__ A text file that includes the coding sequence (CDS) range for all transcripts (see description below). If such an option is not provided, the transcript fasta file is assumed to only include the CDS regions.
* __--work_dir__ (default the parent directory of `scripts`) The working directory where all intermediate and final results will write to.
* __--nproc__ (default 15) Number of threads can be used by ribomap.
* __--adapter__ (default `CTGTAGGCACCATCAAT`) The linker sequence attached to the 5' end of the ribo-seq reads.
* __--nmismatch__ (default 1) Number of mismatches allowed in the read alignments.
* __--softClipping__ (default `true`) Whether reads are allowed to be soft-clipped by STAR when aligning to the transcriptome.
* __--min_fplen__ (default 27) Minimun read length to keep for downstream analysis.
* __--max_fplen__ (default 33) Maximum riboseq read length to keep for downstream analysis.
* __--offset__ (default 12) Offset location in a read that the ribosome P-site maps to, or a text file name that defines the P-site offset based on read length (see description below).
* __--rnaUnstranded__ (default `false`) Whether the RNA-seq protocol is stranded. If the RNA-seq protocol is unstranded, the `librarytype` to run Salifish is set to `-l U`; otherwise the `librarytype` is set to `-l SF`, and alignments with the RC flag set in the RNA-seq data are discarded.
* __--tabd_cutoff__ (default 0) Transcript abundance threshold to be considered expressed.
* __--useSecondary__ (default `true`) Whether multi-mapping alignments are used when assigning footprints to candidate loci.
* __--star_idx_dir__ (default `$work_dir/StarIndex/`) Directory to store Star index.
* __--alignment_dir__ (default  `$work_dir/alignment/`) Directory to store alignment results output by STAR.
* __--sailfish_dir__ (default `$work_dir/sm_quant/`) Directory to store sailfish result.
* __--output_dir__ (default `$work_dir/outputs/`) Directory to store ribomap's outputs.
* __--force__ Force ribomap to regenerate all intermediate steps.

One example of using the shell script:
~~~~~~
    ./run_ribomap.sh \
    --rnaseq_fq rnaseq.fq.gz \
    --riboseq_fq riboseq.fq.gz \
    --contaminant_fa contaminant.fa \
    --transcript_fa transcript.fa \
    --cds_range cds_range.txt
~~~~~~

Please connect the parameter flags and the parameters with a white space.

* __CDS range file__ A plain text file that includes the CDS regions of transcriptome. Each line in the file should be in the following format:

       `transcript_id start stop`

`transcript_id` should be consistent with the transcript ID in the fasta file, `start` is the start base of the coding region in the transcript fasta file (zero-based), and `stop` is one base pass the stop position of the coding sequence. 

* __offset file__ A plain text file that describes the read length and the P site offset for that read length. Each line in the file should be in the following format:

  	   `read_length P-site_offset`

Only a proper range of read length should be included, reads with a length not specified in this file will be discarded for downstream analysis.

#### Run Ribomap by providing the transcript abundance estimation file
Ribomap supports transcript abundance estimation files from [*Sailfish*](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html), [*Cufflinks*](http://cufflinks.cbcb.umd.edu/index.html) and [*eXpress*](http://bio.math.berkeley.edu/eXpress/overview.html). Mapping the ribosome footprint can be performed providing any of the three transcript abundance esitmation files listed above.

Under the `bin` directory, run:

      ./riboprof [options]

The list of options are as follows:

* __-r | --ribobam__ Bam file of ribo-seq read mappings to the transcriptome.
* __-m | --mrnabam__ Bam file of RNA-seq read mappings to the transcriptome.
* __-f | --fasta__ Transcriptome reference fasta file.
* __-cds | --cds_range__ CDS range file.
* __-s | --sf__ Transcript abundance estimation produced by Sailfish.
* __-c | --cl__ Transcript abundance estimation produced by Cufflinks.
* __-e | --ep__ Transcript abundance estimation produced by eXpress.
* __threshold | --tabd_cutoff__ Transcript abundance threshold to be considered expressed.
* __-o | --out__ Output file prefix of Ribomap's result.
* __-p | --offset__ Offset location in a read that the ribosome P-site maps to, or the name of a offset file that specifies P-site offset for different read length.
* __-lmin | --min_fplen__ Minimun read length to keep for downstream analysis.
* __-lmax | --max_fplen__ Maximum riboseq read length to keep for downstream analysis.
* __-sec | --useSecondary__ Use multi-mapping alignments when assigning footprints to candidate loci.
* __-rc | --useRC__ Use alignments with the RC flag set in the RNA-seq data.

One example of using the executable:
~~~~~~
    ./riboprof  \
    --mrnabam mRNA.bam --ribobam ribo.bam \
    --fasta transcript.fa --cds_range cds_range.txt \
    --sf quant.sf  --tabd_cutoff 0 \
    --offset 12  --min_fplen 27 --max_fplen 33 \
    --out ../outputs/ribomap \
    --useSecondary
~~~~~~
Please connect the parameter flags and the parameters with a space.

Ribomap output files
------
Ribomap produces five output files:
#### _XXX.base_
The sub-codon resolution, nucleotide-level ribosome profiles including the UTR regions. Only transcripts with a non-zero total ribosome count are reported. Each entry of a specific transcript looks like this:
~~~~~~
	refID: 0
	tid: YAL001C
	ribo profile: 0 0 0 74 68 ...
	mRNA profile: 31 35 50 73 87 96 104 ...
	normalized ribo profile: 0 0 0 1.0137 0.781609 0.0208333 0.125 ...
~~~~~~  
* __refID__ The transcript fai index in the transcriptome fasta file.
* __tid__ Transcript header name in the transcriptome fasta file.
* __ribo profile__ Nucleotide level ribosome profile including the UTR regions. Each number in the vector is the number of ribosome footprints whose P-sites are covering the corresponding base.
* __mRNA profile__ RNA-seq profile vector of the transcript. Each number in the vector is the read coverage count that are esimated on the corresponding base.
* __normalized ribo profile__ The ribosome profile vector after bias correction. Each number in the vetor is the ratio between the ribo profile count and the mRNA profile count.

#### _XXX.codon_
The in-frame ribosome profiles within the CDS of each transcript.The file format is the same as the the _XXX.base_ file.

#### _XXX.stats_
The summarized statistics for each transcripts. Each entro of a specific transcript looks like this: 
~~~~~
	refID: 0
	tid: YAL001C
	rabd: 3959
	tabd: 0.000209384
	te: 1.89078e+07
~~~~~~
* __refID__ The transcript fai index in the transcriptome fasta file.
* __tid__ Transcript header name in the transcriptome fasta file.
* __rabd__ Ribosome loads, which is the total number of ribosome reads that are esimated from this trascript.
* __tabd__ Relative transcript abundance from Sailfish’s result.
* __te__ Relative translational efficiency, which is the ratio between __rabd__ and __tabd__.


#### _XXX_abundant.list_ 
A list of transcripts whose total ribosome abundance is more than expected given the transcript abundance. 
There is one transcript record per row. The columns are defined as follows:

| Column number | Description |
|---------------|-------------|
| 1 | transcript header name | 
| 2 | relative transcript abundance |
| 3 | total ribosome footprint count |
| 4 | pencentile ranking of the transcript abundance |
| 5 | percentile ranking of the total ribosome footprint count |
| 6 | difference between the transcript abundance rank and the total ribosome footprint count rank

#### _XXX_scarce.list_
A list of transcripts whose total ribosome abundance is less than expected given the transcript abundance.
The file format is the same as _XXX_abundant.list_.

Test case
------
### Run test case
Under the `scripts` directory, run:

      ./hela_ribo_analysis.sh

`hela_ribo_analysis.sh` automatically downloads the transcriptome fasta, gtf, ncRNA fasta, tRNA fasta to the directory `$work_dir/ref/`, and a RNA-seq data and a riboseq data to the directory `$work_dir/fasta/`. The transcriptome fasta file is preprocessed with a _python_ script `filter_gencode_transcript.py` to excludes the following transcripts:

1. transcripts without verified start codon
2. transcripts with stop codon in the middle of the CDS region
3. transcripts with duplicated sequences
4. peptide sequence length less than 3 when the start and stop codons are not included

The CDS range file `gencode.v18.pc_transcripts_cds.txt` is also automatically generated by the same script from parsing the gencode fasta header lines.
A contaminant fasta file `human_contaminant.fa` is built by a _python_ script `build_contaminant.py`. [Biopython](http://biopython.org/wiki/Main_Page) is needed to run the python scripts. 

`hela_ribo_analysis.sh` also sets up the options for `run_ribomap.sh` and automatically calls it.

### Test case data sets
* __RNA-seq__ [GSM546921](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz)
* __riboseq__ [GSM546920](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz)
* __human transcriptome reference fasta__ [gencode.v18.pc_transcripts.fa](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz)
* __human transcriptome annotation gtf__ [gencode.v18.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz)
* __human non-coding rna fasta__ [Homo_sapiens.GRCh38.ncrna.fa] (ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
* __tRNA fasta__ [eukaryotic-tRNAs.fa] (http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz)
