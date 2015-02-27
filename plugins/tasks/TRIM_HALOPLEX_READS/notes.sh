


#!/bin/sh

#
. ${JOB_DIR}/constants.sh

READ_FILES_LIST=""

function plugin_task {
     set -xv
     cd ${TMPDIR}
TRIM_JAR="trimmomatic-0.32.jar"
NUM_THREADS=10

cat >samples.txt <<EOT
1-P1-Re-1
2-P1-Do-1
3-P2-Re-2
4-P2-Do-2
5-P3-Re-3
6-P3-Do-3
7-P4-Re-4
8-P4-Do-4
9-P5-Re-5
10-P5-Do-5
16-P8-Do-8
EOT

cat >samples.txt <<EOT
1-P1-Re-1
2-P1-Do-1
3-P2-Re-2
4-P2-Do-2
5-P3-Re-3
6-P3-Do-3
7-P4-Re-4
8-P4-Do-4
9-P5-Re-5
10-P5-Do-5
11-P6-Re-6
12-P6-Do-6
13-P7-Re-7
14-P7-Do-7
15-P8-Re-8
16-P8-Do-8
17-P9-Re-9
18-P9-Do-9
19-P10-Re-10
20-P10-Do-10
21-P11-Re-11
22-P11-Do-11
23-P12-Re-12
24-P12-Do-12
25-P13-Re-13
26-P13-Do-13
27-P14-Re-14
28-P14-Do-14
29-P15-Re-15
30-P15-Do-15
31-P16-Re-16
32-P16-Do-16
33-P17-Re-17
34-P17-Do-17
35-P18-Re-18
36-P18-Do-18
37-P19-Re-19
38-P19-Do-19
39-P20-Re-20
40-P20-Do-20
41-P21-Re-21
42-P21-Do-21
43-P22-Re-22
44-P22-Do-22
45-P23-Re-23
46-P23-Do-23
47-P24-Re-24
48-P24-Do-24
EOT

mkdir -p cr
mkdir -p trims
rm cr/*
rm trims/*
export TRIM_JAR="trimmomatic-0.32.jar"
cat >script-trim-parallel2.sh <<EOT
SAMPLE="\$1"
PREFIX="\$RANDOM"
find Sample_\${SAMPLE} -name \${SAMPLE}\*_R1\*.fastq.gz >\${PREFIX}_reads.txt;
find Sample_\${SAMPLE} -name \${SAMPLE}\*_R2\*.fastq.gz >\${PREFIX}_pairs.txt;
cat \${PREFIX}_reads.txt;

 parallel -X  --xapply  java -Xmx1g -jar \${TRIM_JAR} PE -threads 1 -phred33 {1} {2} trims/TRIMMED_PAIRED_{1/} trims/TRIMMED_unpaired_{1/} trims/TRIMMED_PAIRED_{2/} trims/TRIMMED_unpaired_{2/} \
ILLUMINACLIP:TruSeq3-PE-revised-Haloplexadapter.fa:3:35:7:5:true \
MAXINFO:40:0.6 LEADING:3 \
TRAILING:3 SLIDINGWINDOW:3:15 HEADCROP:3 CROP:86 MINLEN:40 \
 :::: \${PREFIX}_reads.txt :::: \${PREFIX}_pairs.txt

if [ -e trims/TRIMMED_PAIRED_${SAMPLE}\*_R1\*.fastq.gz -a -e trims/TRIMMED_PAIRED_${SAMPLE}\*_R2\*.fastq.gz ]; then

 java -Xmx1g -jar goby.jar --mode fasta-to-compact --quality-encoding Sanger  --pair-indicator _R1,_R2 --paired-end \
                         trims/TRIMMED_PAIRED_\${SAMPLE}\*_R1\*.fastq.gz \
                         trims/TRIMMED_PAIRED_\${SAMPLE}\*_R2\*.fastq.gz \
                         --concat \
                         -o cr/\${SAMPLE}.compact-reads
fi
rm \${PREFIX}_*.txt
EOT
chmod +x script-trim-parallel2.sh
NUM_THREADS=10
cat samples.txt | parallel --progress  -j ${NUM_THREADS} ./script-trim-parallel2.sh {1}

mkdir -p cr

cat >script-fq-to-cr-parallel.sh <<EOT
SAMPLE="\$1"
PREFIX="\$RANDOM"

find . -name TRIMMED_PAIRED_\${SAMPLE}\*_R1\*.fastq.gz >\${PREFIX}_reads.txt;
find . -name TRIMMED_PAIRED_\${SAMPLE}\*_R2\*.fastq.gz >\${PREFIX}_pairs.txt;

parallel --no-notice -X --xapply echo goby 1g fasta-to-compact --quality-encoding Sanger  --pair-indicator _R1,_R2 --paired-end {1} {2} -o cr/\${SAMPLE}_\${RANDOM}.compact-reads\
    :::: \${PREFIX}_reads.txt :::: \${PREFIX}_pairs.txt

rm \${PREFIX}_*.txt
EOT
chmod +x script-fq-to-cr-parallel.sh

cat samples.txt | parallel --progress --no-notice -j ${NUM_THREADS} ./script-fq-to-cr-parallel.sh {1}



cat >script-fq-to-cr2-parallel.sh <<EOT
SAMPLE="\$1"

goby 1g fasta-to-compact --quality-encoding Sanger  --pair-indicator _R1,_R2 --paired-end \
                         TRIMMED_PAIRED_\${SAMPLE}\*_R1\*.fastq.gz \
                         TRIMMED_PAIRED_\${SAMPLE}\*_R2\*.fastq.gz \
                         --concat \
                         -o cr/\${SAMPLE}.compact-reads

EOT
chmod +x script-fq-to-cr2-parallel.sh

cat samples.txt | parallel --progress --no-notice -j ${NUM_THREADS} ./script-fq-to-cr2-parallel.sh {1}

READS="LM11_AAACATCG_L001_R1_004.fastq.gz"
java -jar trimmomatic-0.32.jar   \
 \
LM11_AAACATCG_L001_R2_004.fastq.gz \
TRIMMED2-Paired-LM11_AAACATCG_L001_R1_004.fastq.gz  \
Trimmed2-Unpaired-LM11_AAACATCG_L001_R1_004.fastq.gz  \
TRIMMED2-Paired-LM11_AAACATCG_L001_R2_004.fastq.gz   \
Trimmed2-Unpaired-LM11_AAACATCG_L001_R2_004.fastq.gz \
 ILLUMINACLIP:TruSeq3-PE-revised-Haloplexadapter.fa:3:35:7:5:true \
MAXINFO:40:0.6 LEADING:3 \
TRAILING:3 SLIDINGWINDOW:3:15 HEADCROP:3 CROP:86 MINLEN:40

     ALL_REGISTERED_TAGS="${LOG_REGISTERED_TAGS} ${VCF_REGISTERED_TAGS}"
     echo "The following tags were registered by this plugin: ${ALL_REGISTERED_TAGS}"
     exit $STATUS
}

    gobyweb@spanky:~/GOBYWEB_SGE_JOBS-dev/fac2003

    ls -1 ?/*/oge*.sh| xargs -n 1 -I{} qsub -terse -l excl=false,h_vmem=11G,virtual_free=11G -r y -v STATE=task {}



java -Xmx1g -jar ${TRIM_JAR} PE -threads 1 -phred33 {1} {2} TRIMMED_PAIRED_{1/} TRIMMED_unpaired_{1/} TRIMMED_PAIRED_{2/} TRIMMED_unpaired_{2/} \
ILLUMINACLIP:TruSeq3-PE-revised-Haloplexadapter.fa:3:35:7:5:true \
MAXINFO:40:0.6 LEADING:3 \
TRAILING:3 SLIDINGWINDOW:3:15 HEADCROP:3 CROP:86 MINLEN:40 \
${PREFIX}_reads.txt ${PREFIX}_pairs.txt

goby 1g fasta-to-compact --quality-encoding Sanger  --pair-indicator _R1,_R2 --paired-end \
                         TRIMMED_PAIRED_\${SAMPLE}\*_R1\*.fastq.gz \
                         TRIMMED_PAIRED_\${SAMPLE}\*_R2\*.fastq.gz \
                         --concat \
                         -o cr/\${SAMPLE}.compact-reads