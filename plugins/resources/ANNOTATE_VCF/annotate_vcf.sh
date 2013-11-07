
function dieUponError {
    RETURN_STATUS=$?
    if [ ! $RETURN_STATUS -eq 0 ]; then
        ls -lat
        # remove the file so that fdr adjustment can proceed (it would fail on an empty file)
        rm $2
        ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "VCF annotate part, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed" --index ${CURRENT_PART} --job-type job-part
    fi

}

function compressWhenNeeded {
  output=$1
  outputNoGzExtension="${output%.gz}"
  if [ ! "${output}" = "${outputNoGzExtension}" ]; then
      ${RESOURCES_TABIX_BGZIP_EXEC_PATH} ${outputNoGzExtension}
  fi
}

function annotate_vep {

    doAnnotate=$1
    input=$2
    output=$3
    doKeepOnlyNonSynonymous=$4
    outputNoGzExtension="${output%.gz}"
    outputTmpVcf="${output%.vcf.gz}-tmp.vcf"
    outputTmpTsv="${output%.vcf.gz}-tmp.tsv"

    inputCopy=`basename "${input%.vcf}-copy.vcf"`


    if [ "${doAnnotate}" == "true" ]; then

        org=`echo ${ORGANISM} | tr '[:upper:]' '[:lower:]'`
        # Disable  --allow_non_variant on April 15 2013 because it produced invalid VCF for sites without annotated variants
        # Instead, we retrieve annotations from vep, extract these annotations to a TSV format, and annotate the original
        # VCF with the new column from the TSV.
        ${RESOURCES_VARIANT_EFFECT_PREDICTOR_SCRIPT} --format vcf -i ${input} -o annotatedInput.vcf --species ${org} \
           --force_overwrite --host useastdb.ensembl.org --vcf
        dieUponError
        export PERL5LIB=${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/lib/perl5/site_perl:${PERL5LIB}
        ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-sort annotatedInput.vcf > sorted.vcf
        ${RESOURCES_TABIX_BGZIP_EXEC_PATH} sorted.vcf
        ${RESOURCES_TABIX_EXEC_PATH} -p vcf sorted.vcf.gz
        # VCF-query refuses to print the same column twice, we need to add this column manually with awk..
        ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-query sorted.vcf.gz -f '%CHROM\t%POS\t%INFO/CSQ\n' >output-dumb.tsv
        awk '{print $1"\t"$2"\t"($2+1)"\t"$3}' output-dumb.tsv >output.tsv;
        #cp output.tsv ${JOB_DIR}/
        ${RESOURCES_TABIX_BGZIP_EXEC_PATH} output.tsv
        ${RESOURCES_TABIX_EXEC_PATH} -p vcf output.tsv.gz

        cp ${input} input.vcf
        ${RESOURCES_TABIX_BGZIP_EXEC_PATH} input.vcf
        ${RESOURCES_TABIX_EXEC_PATH} -p vcf input.vcf.gz

        ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-sort ${input} > raw-input-sorted.vcf
        cat >${TMPDIR}/attributes.lst <<EOF
key=INFO,ID=VariantEffectPrediction,Number=1,Type=String,Description="Variant Effect Predictions"
EOF


         cat raw-input-sorted.vcf | ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-annotate  \
                      -a output.tsv.gz \
                      -d ${TMPDIR}/attributes.lst \
                      -c CHROM,FROM,TO,INFO/VariantEffectPrediction >output-with-vep-info.vcf
         dieUponError
         cp output-with-vep-info.vcf ${JOB_DIR}/
cat >${TMPDIR}/nonSynomymousFilter.pl   <<EOF

# Filter all variants that do not change the coding sequence according to the Variant Effect Prediction Info column.
{
  tag      => 'INFO/VariantEffectPrediction',
  name     => 'NonSynonymousVEP',
  desc     => 'Keep only variants if they are predicted to change a protein sequence',
  apply_to => 'SNPs',
  test     => sub {
                         local \$_ = \$MATCH;
                         if ( /initiator_codon_variant/ ) {return \$PASS; }
                         if ( /inframe_insertion/ ) {return \$PASS; }
                         if ( /inframe_deletion/ ) {return \$PASS; }
                         if ( /missense_variant/ ) {return \$PASS; }
                         if ( /frameshift_variant/ ) {return \$PASS; }
                         if ( /stop_gained/ ) {return \$PASS; }

                         return \$FAIL;

  }
}
EOF
         if [ "${doKeepOnlyNonSynonymous}" == "true" ]; then
            cat output-with-vep-info.vcf | ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-annotate \
                      -f nonSynomymousFilter.pl \
                      >${outputNoGzExtension}
            dieUponError
         else
           cp output-with-vep-info.vcf ${outputNoGzExtension}
         fi

        compressWhenNeeded ${output}
        dieUponError

    else

        cp ${input} ${outputNoGzExtension}
        dieUponError
        compressWhenNeeded ${output}
        dieUponError
    fi

}


function annotate_ensembl_genes {

    doAnnotate=$1
    input=$2
    output=$3
    outputNoGzExtension="${output%.gz}"

    if [ "${doAnnotate}" == "true" ]; then
       ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
       BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
       ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
       ANNOTATION_PATH=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_ANNOTATIONS_ANNOTATIONS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

           cat >${TMPDIR}/attributes.lst <<EOF
key=INFO,ID=GENE,Number=1,Type=String,Description="Ensembl gene identifier"
key=INFO,ID=GENE_NAME,Number=1,Type=String,Description="Gene name"
EOF

     ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-sort ${input} \
       | ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-annotate -a ${ANNOTATION_PATH}/ref-start-end-gene-hgnc-sorted.tsv.gz \
                      -d ${TMPDIR}/attributes.lst -c CHROM,FROM,TO,INFO/GENE,INFO/GENE_NAME >${outputNoGzExtension}

     dieUponError

     compressWhenNeeded ${output}
     dieUponError


    else

        cp ${input} ${outputNoGzExtension}
        dieUponError
        compressWhenNeeded ${output}
        dieUponError
    fi

}

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

     ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-sort ${input} \
       | ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-annotate -a ${TMPDIR}/vep-results.tsv.gz \
                      -d ${TMPDIR}/vep.lst -c CHROM,FROM,TO,INFO/GENE,INFO/GENE_NAME,INFO/RS,INFO/Effect  \
       |${RESOURCES_TABIX_BGZIP_EXEC_PATH} -c > ${output}

}