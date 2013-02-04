#!/bin/sh

function dieUponError {
    RETURN_STATUS=$?
    if [ ! $RETURN_STATUS -eq 0 ]; then
        ls -lat
        # remove the file so that fdr adjustment can proceed (it would fail on an empty file)
        rm $2
        %QUEUE_WRITER% --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "VCF annotate part, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed" --index ${CURRENT_PART} --job-type job-part
    fi

}

function annotate {

    input=$1
    output=$2
    outputNoGzExtension="${output%.gz}"
    org=`echo ${ORGANISM} | tr '[:upper:]' '[:lower:]'`
    # Retrieve annotations from vep and rewrite the VCF:
    ${RESOURCES_VARIANT_EFFECT_PREDICTOR_SCRIPT} --format vcf -i ${input} -o ${outputNoGzExtension} --species ${org} \
       --force_overwrite --host useastdb.ensembl.org --vcf
    dieUponError

    # make extension
    ${RESOURCES_TABIX_BGZIP_EXEC_PATH} ${outputNoGzExtension}
    dieUponError

}

annotate $1 $2

function dummy_ignore_this {
  ${RESOURCES_VARIANT_EFFECT_PREDICTOR_SCRIPT} --format vcf -i ${input} -o ${TMPDIR}/vep-results.tsv --species ${org} \
        --force_overwrite --host useastdb.ensembl.org


    # index with bgzip and tabix:
    ${RESOURCES_TABIX_BGZIP_EXEC_PATH} vep-results.tsv
    ${RESOURCES_TABIX_EXEC_PATH} vep-results.tsv.gz

    #Combine to output with vcf-annotate:
         cat >${TMPDIR}/vep.lst <<EOF
key=INFO,ID=GENE,Number=1,Type=String,Description="Ensembl gene identifier"
key=INFO,ID=GENE_NAME,Number=1,Type=String,Description="Gene name"
key=INFO,ID=RS,Number=1,Type=String,Description="rs ID"
key=INFO,ID=Effect,Number=1,Type=String,Description="Effect of variation on transcript, predicted with VEP"
EOF

     ${VCFTOOLS_BIN}/vcf-sort ${input} \
       | ${VCFTOOLS_BIN}/vcf-annotate -a ${TMPDIR}/vep-results.tsv.gz \
                      -d ${TMPDIR}/vep.lst -c CHROM,FROM,TO,INFO/GENE,INFO/GENE_NAME,INFO/RS,INFO/Effect  \
       |${BGZIP_EXEC_PATH} -c > ${output}

}