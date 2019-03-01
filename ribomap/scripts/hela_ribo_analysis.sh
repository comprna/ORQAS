#!/bin/bash
src_dir=`dirname $0`
ribomap_dir=${src_dir}/../
work_dir=${ribomap_dir}
fasta_dir=${work_dir}data/fasta/
ref_dir=${work_dir}ref/
get_data=true
#=============================
# urls
#=============================
riboseq_url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz
rnaseq_url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz
gtf_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz
tfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz
pfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_translations.fa.gz
gtf=${ref_dir}${gtf_url##*/}
gtf=${gtf%.gz}
tfa=${ref_dir}${tfa_url##*/}
tfa=${tfa%.gz}
pfa=${ref_dir}${pfa_url##*/}
pfa=${pfa%.gz}
if [ "${get_data}" = true ]; then
    #=============================
    # step 1: download sra & refs
    #=============================
    echo "downloading Hela cell reads..."
    wget -P ${fasta_dir} -N ${rnaseq_url}
    wget -P ${fasta_dir} -N ${riboseq_url}
    echo "downloading transcriptome reference data..."
    wget -P ${ref_dir} -N ${gtf_url}
    wget -P ${ref_dir} -N ${tfa_url}
    wget -P ${ref_dir} -N ${pfa_url}
    echo "unzipping data..."
    gunzip -f ${ref_dir}*.gz
    #=============================
    # step 2: process reference
    # build cds range file
    #=============================
    echo "filtering transcripts...."
    python ${src_dir}/filter_gencode_transcript.py ${gtf} ${tfa} ${pfa}
    #====================================
    # step 3: build contaminant sequence
    #====================================
    nc_url=ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
    trna_url=http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
    echo "downloading ncRNA from Ensembl..."
    wget -P ${ref_dir} -N ${nc_url}
    echo "downloading tRNA from gtrnadb..."
    wget -P ${ref_dir} -N ${trna_url}
    gunzip -f ${ref_dir}*.gz
    echo "merging rRNA and tRNA..."
    rrna_fa=${ref_dir}${nc_url##*/}
    rrna_fa=${rrna_fa%.gz}
    trna_fa=${ref_dir}${trna_url##*/}
    trna_fa=${trna_fa%.gz}
    python ${src_dir}/build_contaminant.py ${rrna_fa} ${trna_fa} Homo_sapiens ${ref_dir}human_contaminant.fa
    echo "done preparing data for ribomap"
fi
#=============================
# step 4: run ribomap
#=============================
rnaseq_fq=${fasta_dir}${rnaseq_url##*/}
riboseq_fq=${fasta_dir}${riboseq_url##*/}
transcript_fa=${tfa%.*}_filter.fa
contaminant_fa=${ref_dir}human_contaminant.fa
cds_range=${tfa%.*}_cds.txt
adapter=TCGTATGCCGTCTTCTGCTTG
min_fplen=25
max_fplen=36
offset=${src_dir}/offset.txt
ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --adapter ${adapter} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --offset ${offset}"
# ribomap
${ribo_cmd} --output_dir ${work_dir}outputs #--force true
# star prime
# ${ribo_cmd} --output_dir ${work_dir}star_prime --tabd_cutoff -1 --useSecondary false
