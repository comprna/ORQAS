#############################################
### ORQAS PIPELINE                        ###
### Computational RNA Biology Group       ###
### Marina Reixachs  2019                 ###
#############################################




#EXAMPLE COMMAND: 
#sh run_ORQAS.sh --rnaseq_fq rnaseq_fq --riboseq_fq riboseq_fq --cds_annotation cds.fa --txt_cds cds_annotation.ENSTtoCDS.txt --salmon_idx index_dir --salmon_strand U --psites psites.txt --wd orqas_output_directory

#to provide
while [[ $# > 1 ]]
do
    key="$1"
    shift
    case $key in
	--rnaseq_fq)
	    rnaseq_fq="$1"
	    shift
	    ;;
	--riboseq_fq)
	    riboseq_fq="$1"
	    shift
	    ;;
	--orqas_dir)
	    orqas_dir="$1"
	    shift
	    ;;
	--cds_annotation)
	    cds_annotation="$1"
	    shift
	    ;;
	--txt_cds)
	    txt_cds="$1"
	    shift
	    ;;
	--salmon_idx)
	    salmon_idx="$1"
	    shift
	    ;;
	--salmon_strand)
	    salmon_strand="$1"
	    shift
	    ;;
	--psites)
	    psites="$1"
	    shift
	    ;;
	--wd)
	    wd="$1"
	    shift
	    ;;
    esac
done

if [ -z "${orqas_dir}" ]; then
    orqas_dir="." 
fi

if [ -z "${riboseq_fq}" ]; then 
    error_msg "ribo-seq reads not provided!" true
elif [ ! -f ${riboseq_fq} ]; then
    error_msg "ribo-seq file not exist! ${riboseq_fq}" true
elif [ -z "${rnaseq_fq}" ]; then
    error_msg "RNA-seq reads not provided!" true
elif [ ! -f ${rnaseq_fq} ]; then
    error_msg "RNA-seq file not exist! ${rnaseq_fq}" true
elif [ -z "${salmon_idx}" ]; then
    salmon_idx="${wd}/sm_index/"
     ${orqas_dir}/salmon-0.7.2/bin/salmon index -t ${cds_annotation} -i ${salmon_idx} --type quasi -k 31
    error_msg "Salmon index not provided! A new one with default parameters will be generated." true
elif [ ! -d ${salmon_idx} ]; then
    salmon_idx="${wd}/sm_index/"
    ${orqas_dir}/salmon-0.7.2/bin/salmon index -t ${cds_annotation} -i ${salmon_idx} --type quasi -k 31
    error_msg "Salmon index not exists ${salmon_idx}! A new one with default parameters will be generated." true
fi





###################
# 1: RUN SALMON   #
###################

salmon_output="${wd}/sm_quant/"



if [[ -e  ${salmon_output} ]]; then
      echo "Salmon output already exists. Skipping..."
else 
      ${orqas_dir}/salmon-0.7.2/bin/salmon quant -r ${rnaseq_fq} -i ${salmon_idx} --libType ${salmon_strand} -o ${salmon_output}
fi


####################
# 2: RUN RIBOMAP   #
####################

if [[ ${salmon_strand} == *"U"* ]]; then
     rnaseq_strand="true"
else
     rnaseq_strand="false"
fi

	
if [[ -e  "${wd}/outputs/" ]]; then
      echo "Ribomap output already exists. Skipping..."
else 
      ${orqas_dir}/ribomap/scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${cds_annotation} --work_dir "${wd}/" --offset ${psites} --softClipping false --rnaUnstranded ${rnaseq_strand}  --sailfish_dir ${salmon_output}
fi


######################
# 3: CALCULATE OPM   #
######################

name=$(basename ${wd})
out="${wd}/${name}"

if [[ -e "${out}.Abundance.txt" ]]; then
      echo "OPMcalculator output already exists. Skipping..."
else 
	python ${orqas_dir}/ORQAStoolkit/ORQAStoolkit.py OPMcalculator -i "${wd}/" -o "${out}"
fi


######################
# 4: AGGREGATE CDS   #
######################

if [[ -e  "${wd}/${name}.Abundance.aggregateCDS.txt" ]]; then
      echo "aggregateCDS output already exists. Skipping..."
else

     if [ -z "${txt_cds}" ]; then 
          echo "Transcript to cds assignment not provided... calculating it."
          txt_cds="${wd}/${name}"
          python ${orqas_dir}/ORQAStoolkit/ORQAStoolkit.py TXTtoCDS -i ${cds_annotation} -o ${txt_cds}
     fi

     python ${orqas_dir}/ORQAStoolkit/ORQAStoolkit.py aggregateCDS -i "${wd}/${name}.Abundance.txt" -o "${wd}/${name}.Abundance.aggregateCDS.txt" -c "${txt_cds}.ENSTtoCDS.txt" 
fi

echo "DONE! Outputs are in ${wd}"
echo "To obtain periodicity and uniformity values run validateORF with a single sample or a pool of samples from the same condition"


