#!/bin/env groovy

@Grab(group = 'org.apache.httpcomponents', module = 'httpclient', version = '4.1.1')

import com.martiansoftware.jsap.JSAPResult
import org.apache.http.HttpEntity
import org.apache.http.HttpResponse
import org.apache.http.NameValuePair
import org.apache.http.client.HttpClient
import org.apache.http.client.entity.UrlEncodedFormEntity
import org.apache.http.client.methods.HttpPost
import org.apache.http.impl.client.DefaultHttpClient
import org.apache.http.message.BasicNameValuePair
import org.campagnelab.groovySupport.CurlSupport
import org.campagnelab.groovySupport.ExecAndRemote

/**
 * Groovy utility to download tab delimited data from Biomart.
 * -e Specifies the export to run, can 'all' or a comma separated list
 *    from exon-annotations,var-annotations,ref-start-end-gene. Default is 'exon-annotations'
 * -f The output folder to write files to. Default is '.'.
 * -r The organism to retrieve from homo_sapiens,mus_musculus,canis_familiaris,m_tuberculosis_h37rv,danio_rerio,rabbit
 * -d Each export above in -e has a "default" dataset. If you to override the dataset,
 *    specify it here.
 * -v The virtual schema name. Default is 'default'.
 * -u The server prefix. Defaul is 'uswest'.
 * -h Help
 *
 * On the tabix files, it can be useful to run a sanity check such as
 *     tabix ref-start-end-gene-sorted.tsv.gz 11
 * to check that things went well (chromosome 11 should be retrieved, in this case). If no output
 * is returned the file is bad.
 *
 * For turberculosis exon-annotations, the dataset needs to be changed. Example:
 *    ./Biomart.groovy -r m_tuberculosis_h37rv -d 30_gene -f test-turb -e exon-annotations
 *
 */
public class Biomart {

    final String XML_CONFIG =
        """
<jsap>
    <parameters>
        <flaggedOption>
            <id>exports</id>
            <longFlag>exports</longFlag>
            <shortFlag>e</shortFlag>
            <required>false</required>
            <defaults>
                <string>%EXPORTS_DEFAULT%</string>
            </defaults>
            <help>Which exports to run, either 'all' or a comma separated of %EXPORTS_LIST%.</help>
        </flaggedOption>
        <flaggedOption>
            <id>output-folder</id>
            <longFlag>output-folder</longFlag>
            <shortFlag>f</shortFlag>
            <required>false</required>
            <defaults>
                <string>.</string>
            </defaults>
            <help>The folder to write the output file(s) to.</help>
        </flaggedOption>
        <flaggedOption>
            <id>organism</id>
            <longFlag>organism</longFlag>
            <shortFlag>r</shortFlag>
            <required>false</required>
            <defaults>
                <string>%ORGANISMS_DEFAULT%</string>
            </defaults>
            <help>Organism, one of %ORGANISMS_LIST%.</help>
        </flaggedOption>
        <flaggedOption>
            <id>dataset</id>
            <longFlag>dataset</longFlag>
            <shortFlag>d</shortFlag>
            <required>false</required>
            <help>Biomart dataset to query.</help>
        </flaggedOption>
        <flaggedOption>
            <id>virtual-schema-name</id>
            <longFlag>virtual-schema-name</longFlag>
            <shortFlag>v</shortFlag>
            <required>false</required>
            <defaults>
                <string>default</string>
            </defaults>
            <help>The Biomart/Martservice virtual-schema-name to use.</help>
        </flaggedOption>
        <flaggedOption>
            <id>url-prefix</id>
            <longFlag>url-prefix</longFlag>
            <shortFlag>u</shortFlag>
            <required>false</required>
            <defaults>
                <string>uswest</string>
            </defaults>
            <help>URL prefix for archive (such as uswest or useast for current data or such as jul2009 for historical data).</help>
        </flaggedOption>
        <flaggedOption>
            <id>chromosome-list-file</id>
            <longFlag>chromosome-list-file</longFlag>
            <shortFlag>c</shortFlag>
            <required>false</required>
            <help>A text file listing chromosomes, one per line.</help>
            <stringParser>
                <classname>FileStringParser</classname>
            </stringParser>
        </flaggedOption>
        <flaggedOption>
            <id>tabix-executable</id>
            <longFlag>tabix-executable</longFlag>
            <required>true</required>
            <help>The path to the tabix executable.</help>
            <stringParser>
                <classname>FileStringParser</classname>
                <properties>
                    <property>
                        <name>mustBeDirectory</name>
                        <value>false</value>
                    </property>
                    <property>
                        <name>mustExist</name>
                        <value>true</value>
                    </property>
                </properties>
            </stringParser>
        </flaggedOption>
        <flaggedOption>
            <id>bgzip-executable</id>
            <longFlag>bgzip-executable</longFlag>
            <required>true</required>
            <help>The path to the bgzip executable.</help>
            <stringParser>
                <classname>FileStringParser</classname>
                <properties>
                    <property>
                        <name>mustBeDirectory</name>
                        <value>false</value>
                    </property>
                    <property>
                        <name>mustExist</name>
                        <value>true</value>
                    </property>
                </properties>
            </stringParser>
    </flaggedOption>
</parameters>
</jsap>
"""

    def env = System.getenv()
    def TABIX_EXEC = null
    def BGZIP_EXEC = null
    def exec = new ExecAndRemote()
    def curl = new CurlSupport()

    boolean deleteUnsorted = true

    boolean downloadWithCurl = false

    boolean progressTicks = true

    boolean urlEncodeQuery = true

    boolean forceNoFilterByChrom = false

    /** This can be useful for re-fetching chromosomes that failed. */
    boolean omitAllHeaders = false

    /** This can be useful for re-fetching chromosomes that failed.*/
    def limitFetchToChromosomesList = []

    int maxRetries = 1

    // "m_tuberculosis_h37rv this will require -d 30_gene
    static organismToDatasetPrefix = [
            "homo_sapiens": "hsapiens",
            "mus_musculus": "mmusculus",
            "canis_familiaris": "cfamiliaris",
            "m_tuberculosis_h37rv": "myc",
            "danio_rerio": "drerio",
            "rattus_norvegicus": "rnorvegicus",
            "xenopus_tropicalis": "xtropicalis",
            "Oryctolagus_cuniculus": "ocuniculus",
            "Caenorhabditis_elegans":"celegans",

    ]

    static exportTypeToDataMap = [
            "exon-annotations": [
                    outputFilename: "exon-annotations.tsv",
                    sort: "-k 1,1 -k 5,5n -k 6,6n",
                    dataset: "gene_ensembl",
                    fields: [
                            "chromosome_name": "Chromosome Name",
                            "strand": "Strand",
                            "ensembl_gene_id": "Ensembl Gene ID",
                            "ensembl_exon_id": "Ensembl Exon ID",
                            "exon_chrom_start": "Exon Chr Start (bp)",
                            "exon_chrom_end": "Exon Chr End (bp)",
                    ]
            ],
            "gene-annotations": [
                                          outputFilename: "gene-annotations.tsv",
                                          sort: "-k 1,1 -k 5,5n -k 6,6n",
                                          dataset: "gene_ensembl",
                                          fields: [
                                                  "chromosome_name": "Chromosome Name",
                                                  "strand": "Strand",
                                                  "ensembl_gene_id": "Ensembl Gene ID",
                                                  "ensembl_gene_id": "Ensembl Gene ID",
                                                  "start_position": "FROM",
                                                  "end_position": "TO",
                                          ]
                                  ],
           "five-prime-annotations": [
                               outputFilename: "five-prime-annotations.tsv",
                               sort: "-k 1,1 -k 5,5n -k 6,6n",
                               dataset: "gene_ensembl",
                               fields: [
                                       "chromosome_name": "Chromosome Name",
                                       "strand": "Strand",
                                       "ensembl_gene_id": "Ensembl Gene ID",
                                       "ensembl_transcript_id": "Ensembl Transcript ID",
                                       "5_utr_start": "FROM",
                                       "5_utr_end": "TO",
                                ]
                       ],
            "three-prime-annotations": [
                                           outputFilename: "three-prime-annotations.tsv",
                                           sort: "-k 1,1 -k 5,5n -k 6,6n",
                                           dataset: "gene_ensembl",
                                           fields: [
                                                   "chromosome_name": "Chromosome Name",
                                                   "strand": "Strand",
                                                   "ensembl_gene_id": "Ensembl Gene ID",
                                                   "ensembl_transcript_id": "Ensembl Transcript ID",
                                                   "3_utr_start": "FROM",
                                                   "3_utr_end": "TO",
                                            ]
                                   ],
            "ref-start-end-gene": [
                    outputFilename: "ref-start-end-gene-sorted.tsv",
                    sort: "-k 1,1 -k 2,2n -k 3,3n",
                    index: "-s 1 -b 2 -e 3",
                    dataset: "gene_ensembl",
                    fields: [
                            "chromosome_name": "CHROM",
                            "start_position": "FROM",
                            "end_position": "TO",
                            "ensembl_gene_id": "GENE",
                    ]
            ],
            "ref-start-end-gene-hgnc": [
                                outputFilename: "ref-start-end-gene-hgnc-sorted.tsv",
                                sort: "-k 1,1 -k 2,2n -k 3,3n",
                                index: "-s 1 -b 2 -e 3",
                                dataset: "gene_ensembl",
                                fields: [
                                        "chromosome_name": "CHROM",
                                        "start_position": "FROM",
                                        "end_position": "TO",
                                        "ensembl_gene_id": "GENE",
                                        "hgnc_symbol": "HGNC",
                                ]
                        ],
            "gene-id-description": [
                    outputFilename: "gene_id_description.tsv",   // note the undescores because we make a SqlLite db
                    // file and dashes in the filename would prevent attaching the database
                    database: true,
                    dataset: "gene_ensembl",
                    fields: [
                            "ensembl_gene_id": "GENE",
                            "external_gene_id": "GENE_ID",  // a gene id/gene name
                            "description": "DESCRIPTION",
                    ]
            ],
            "var-annotations": [
                    outputFilename: "var-annotations-sorted.tsv",
                    sort: "-k 1,1 -k 2,2n -k 3,3n",
                    index: "-s 1 -b 2 -e 2",
                    dataset: "snp",
                    filterByChrom: true,
                    optional: true, /** Indicate that failure to retrieve this dataset should not stop the import.
             Some organisms do not have variation databases.  */
                    fields: [
                            "chr_name": "CHROM",
                            "chrom_start": "POS",
                            "refsnp_id": "RS",
                            "ensembl_gene_stable_id": "GeneID",
                            "consequence_type_tv": "EFFECT",
                    ]
            ],
    ]

    // Usually "default" but for barteria may be something like "bacterial_mart_7" for release 7 of bacteria
    String virtualSchemaName

    // www for current or something else for archives, such as jul2009
    String serverPrefix

    String organism

    String[] exportsList

    String outputFolder

    String commandLineDataset

    File chromListFile

    public static void main(final String[] args) {
        new Biomart().execute(args)
    }

    def execute(String[] args) {
        if (!configure(args)) {

            System.exit(1)
        }
        for (String exportType in exportsList) {
            int retries = 0
            boolean completed = false
            while (!completed && retries < maxRetries) {
                // Sometimes ensembl returns an error, retry 3 times before giving up
                if (fetch(exportType)) {
                    completed = true
                    break
                } else {
                    if (exportTypeToDataMap[exportType].optional) {
                        System.err.println("The previous error will be ignored because the dataset ${exportType} is optional.")

                        //failure to fetch an optional file is not an error.
                        completed = true;
                        break;
                    }

                    retries++
                    System.err.println "FETCH DETAIL # ${retries} FAILED. Will retry up to 3 times."
                }
            }
            if (!completed) {
                System.err.println "Could not fetch data for ${exportType}, tried 3 times. Aborting."
                System.exit(1)
            }
        }
    }
    /**
     *
     * @param exportType
     * @return success: False when failure.
     */
    def fetch(final String exportType) {
        Map config = exportTypeToDataMap[exportType]

        String dataset = "${organismToDatasetPrefix[organism]}_${commandLineDataset ?: config.dataset}"
        Map fields = config.fields
        String outputFilename = "${this.outputFolder}/${config.outputFilename}"
        String sortedOutputFilename = null
        boolean filterByChrom = config.filterByChrom ?: false
        if (forceNoFilterByChrom) {
            filterByChrom = false
        }
        if (config.sort) {
            sortedOutputFilename = outputFilename
            // Try to replace "sorted" with "unsorted"
            outputFilename = outputFilename.replace("sorted", "unsorted")
            if (outputFilename == sortedOutputFilename) {
                // Doesn't contain "sorted", just append "-unsorted"
                outputFilename = sortedOutputFilename + "-unsorted"
            }
        }
        // def url = "http://${serverPrefix}.ensembl.org/biomart/martservice"
        def url = "http://www.biomart.org/biomart/martservice"
        System.err.println "--------------------------------------------"
        System.err.println "exportType=${exportType}"
        System.err.println "virtualSchemaName=${virtualSchemaName}"
        System.err.println "serverPrefix=${serverPrefix}"
        System.err.println "dataset=${dataset}"
        System.err.println "outputFilename=${outputFilename}"
        System.err.println "fields=${fields}"
        System.err.println "url=${url}"

        def chromList = null
        if (filterByChrom) {
            if (chromListFile == null) {
                println "filterByChrome set to true but no chromosome-list-file specified"
                return
            } else {
                System.err.println "chromListFile=${chromListFile}"
            }
            // Build the chrom list from the file
            chromList = []
            chromListFile.eachLine { chrom ->
                if (chrom) {
                    chromList << chrom
                }
            }
            System.err.println "chromList=${chromList}"
        }

        if (chromList == null) {
            if (!fetchFile(url, outputFilename, true, virtualSchemaName, dataset, fields, null)) {
                // failure to fetch an optional file is not an error:
                return !config.optional
            }
        } else {
            int i = 0
            int numChroms = chromList.size()
            def filesToMerge = []
            for (String chrom in chromList) {
                String chromFilename = outputFilename + "." + chrom
                filesToMerge << "'" + chromFilename + "'"
                boolean writeHeader
                if (omitAllHeaders) {
                    writeHeader = false
                } else {
                    writeHeader = (i++ == 0)
                }

                boolean fetchCurrentFile = true
                if (limitFetchToChromosomesList && (!limitFetchToChromosomesList.contains(chrom))) {
                    fetchCurrentFile = false
                }

                if (fetchCurrentFile) {
                    println ">"
                    println "> Fetching chromosome ${chrom}, ${i} of ${numChroms} for ${organism} / ${exportType}"
                    println ">"
                    if (!fetchFile(url, chromFilename, writeHeader, virtualSchemaName, dataset, fields, chrom)) {
                        return !config.optional
                    }
                }
            }
            exec.bashExec "cat ${filesToMerge.join(" ")} > ${outputFilename}"
            if (deleteUnsorted) {
                exec.exec "rm ${filesToMerge.join(" ")}"
            }
        }

        if (config.index) {
            config.compress = true
        }
        if (config.sort) {
            // Retain the header as the first line when sorting
            exec.bashExec "sed '1q' '${outputFilename}' > '${sortedOutputFilename}'"
            exec.bashExec "sed '1d' '${outputFilename}' | sort ${config.sort} >> '${sortedOutputFilename}'"
        }
        if (config.compress) {
            // Index files MUST be compressed
            String uncompressedFilename = "${sortedOutputFilename ?: outputFilename}"
            String compressedFilename = "${uncompressedFilename}.gz"
            File compressedFile = new File(compressedFilename)
            if (compressedFile.exists()) {
                compressedFile.delete()
            }
            exec.exec "${BGZIP_EXEC} '${uncompressedFilename}'"
            sortedOutputFilename = sortedOutputFilename ? "${sortedOutputFilename}.gz" : null
            outputFilename += ".gz"
        }
        if (config.index) {
            String unindexedFilename = "${sortedOutputFilename ?: outputFilename}"
            String indexedFilename = "${unindexedFilename}.tbi"
            File indexedFile = new File(indexedFilename)
            if (indexedFile.exists()) {
                indexedFile.delete()
            }
            exec.exec "${TABIX_EXEC} ${config.index} '${unindexedFilename}'"
        }
        if (config.database) {
            TsvVcfToSqlite importer = new TsvVcfToSqlite();
            importer.addInputFilename(outputFilename)
            importer.tag = "sdhkshd82392wniks"  // install a dummy tag to prevent replacing "_" with ""

            int status = importer.execute()
            if (status != 0) {
                new File(importer.getOutputFilename()).delete();
                System.err.println("Database import failed for ${organism} ${outputFilename} status=${status}");
                System.err.println("Import file was ${outputFilename}, tried to write to ${importer.getOutputFilename()}");
            }
        }
        if (config.sort) {
            // Delete the unsorted version
            if (deleteUnsorted) {
                new File(outputFilename).delete()
                new File(outputFilename - ".gz").delete()
            }
        }

        System.err.println "exportType export to ${sortedOutputFilename ?: outputFilename} completed."
        return true
    }

    /**
     *
     * @param url
     * @param outputFilename
     * @param writeHeader
     * @param virtualSchemaName
     * @param dataset
     * @param fields
     * @param filterChrom
     * @return
     */
    def fetchFile(url, outputFilename, writeHeader, virtualSchemaName, dataset, fields, filterChrom) {
        // Output the header
        PrintStream out = null
        File xmlQueryFile = null
        File outputFile = new File(outputFilename)
        try {
            // Generate the XML for the fetch
            StringBuffer xmlsb = new StringBuffer()
            xmlsb << """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>"""
            xmlsb << """<Query virtualSchemaName="${virtualSchemaName}" formatter="TSV" header="0" """
            xmlsb << """uniqueRows="0" count="" datasetConfigVersion="0.6" >"""

            xmlsb << """<Dataset name="${dataset}" interface = "default" >"""
            if (filterChrom) {
                def filterField = (fields.keySet() as List)[0]
                xmlsb << """<Filter name="${filterField}" value="${filterChrom}"/>"""
            }
            for (String field in fields.keySet()) {
                xmlsb << """<Attribute name="${field}"/>"""
            }
            xmlsb << """</Dataset></Query>"""
            String xml = xmlsb.toString()
            System.err.println("XML=${xml}")

            if (downloadWithCurl) {
                def now = new Date()
                xmlQueryFile = new File("TEMP-MARTSERVICE-QUERY-" + String.format('%tF', now) + "-" + String.format('%tH-%<tM-%<tS.%<tL', now) + ".xml")
                xmlQueryFile.write "query=${urlEncodeQuery ? URLEncoder.encode(xml) : xml}"

                // Write the header
                if (writeHeader) {
                    File headerFile = new File("${outputFilename}.header")
                    File bodyFile = new File("${outputFilename}.body")
                    PrintWriter headerWriter = new PrintWriter(headerFile.newWriter())
                    headerWriter.println fields.values().join("\t")
                    headerWriter.close()
                    int exitCode = exec.exec("curl -o ${bodyFile} --tr-encoding -d @${xmlQueryFile} ${url}")
                    String exitStatus = curl.exitStatus(exitCode)
                    println "File fetch status: ${exitStatus}"
                    if (exitCode != 0) {

                        return false
                    }
                    if (!bodyFile.exists()) {
                        // If martservice had no data for the chromosome, curl won't write a file. Make the file.
                        exec.exec "touch ${bodyFile}"
                    }
                    exec.bashExec "cat '${headerFile}' '${bodyFile}' > '${outputFilename}'"
                    if (deleteUnsorted) {
                        bodyFile.delete()
                        headerFile.delete()
                    }
                } else {
                    int exitCode = exec.exec("curl -o ${outputFilename} --tr-encoding -d @${xmlQueryFile} ${url}")
                    String exitStatus = curl.exitStatus(exitCode)
                    println "File fetch status: ${exitStatus}"
                    if (exitCode != 0) {
                        return false
                    }
                    if (!outputFile.exists()) {
                        // If martservice had no data for the chromosome, curl won't write a file. Make the file.
                        exec.exec "touch ${outputFilename}"
                    }
                }
            } else {
                out = new PrintStream(outputFile)

                // Output the header
                if (writeHeader) {
                    out.println fields.values().join("\t")
                }

                // Post the data
                List<NameValuePair> formparams = new ArrayList<NameValuePair>();
                formparams.add(new BasicNameValuePair("query", xml));
                UrlEncodedFormEntity formdata = new UrlEncodedFormEntity(formparams, "UTF-8");
                HttpClient httpclient = new DefaultHttpClient();
                HttpPost httpPost = new HttpPost(url);
                httpPost.setEntity(formdata);

                // Retrieve the response
                HttpResponse response = httpclient.execute(httpPost);
                HttpEntity responseEntity = response.getEntity();
                Date lastTime=new Date()
                if (responseEntity != null) {
                    InputStream instream = responseEntity.getContent();
                    int l
                    long total = 0
                    def bufferSize = 10 * 1024 * 1024
                    byte[] tmp = new byte[bufferSize]
                    while ((l = instream.read(tmp)) != -1) {
                        out.print(new String(tmp, 0, l))
                        total += l
                        if (progressTicks) {
                            //print "."
                            Date currentTime = new Date()
                            long elapsedSecs = (currentTime.getTime() - lastTime.getTime()) / 1000;
                            float bytesPerSec = ((float)l) / (elapsedSecs+1);
                            float kbPerSec = ((float)l) / 1024f / (elapsedSecs+1);
                            float mbPerSec = ((float)l) / 1024f / 1024f / (elapsedSecs+1);
                            printf("Downloaded %d so far at %f bytes/sec %f kb/sec %f MB/sec\r", total,
                                    bytesPerSec,kbPerSec, mbPerSec)
                            lastTime = currentTime
                        }
                    }
                    println "Received ${total} at ${new Date().toString()}"
                }
            }
        } catch (Exception e){
            e.printStackTrace()
            System.exit(127);
        } finally{
            if (out != null) {
                out.close()
            }
            if (xmlQueryFile != null && deleteUnsorted) {
                xmlQueryFile.delete()
            }
        }
        return validateFile(outputFile, filterChrom)
    }

    /**
     *
     * @param file
     * @param filterChrom
     * @return False upon error, True upon success
     */
    def validateFile(File file, filterChrom) {
        if (file.length() == 0) {
            if (filterChrom == null) {
                System.err.println "File length was 0"
                return false
            } else {
                println "Chromosome ${filterChrom} contained no data"
            }
        }
        if (file.length() < 2048) {
            // Very small file, let's check the contents for an error
            String fileContents = file.text
            if (fileContents.contains("<title>302 Found</title>")) {
                System.err.println "File contained '<title>302 Found</title>'. Some error in url or post."
                System.err.println ""
                System.err.println "File Contents:"
                System.err.println "-------------------"
                System.err.println fileContents
                return false
            }
            if (fileContents.contains("Query ERROR")) {
                System.err.println "File contained 'Query ERROR'. Some error in url or post."
                System.err.println ""
                System.err.println "File Contents:"
                System.err.println "-------------------"
                System.err.println fileContents
                return false
            }
        }
        return true
    }
    /**
     *
     * @param args
     * @return False upon error, True upon success
     */
    def boolean configure(final String[] args) {
        def jsapConfig = new org.campagnelab.groovySupport.JsapSupport()
                .setArgs(args)
                .setXmlConfig(XML_CONFIG)
                .setScriptName("Biomart.groovy")
                .setHelpValues(
                ["%EXPORTS_LIST%": exportTypeToDataMap.keySet().join(", "),
                        "%EXPORTS_DEFAULT%": exportTypeToDataMap.keySet().asList()[0],
                        "%ORGANISMS_LIST%": organismToDatasetPrefix.keySet().join(", "),
                        "%ORGANISMS_DEFAULT%": organismToDatasetPrefix.keySet().asList()[0],
                ])
        JSAPResult jsap = jsapConfig.parse()

        TABIX_EXEC = jsap.getFile("tabix-executable").getAbsolutePath()
        BGZIP_EXEC = jsap.getFile("bgzip-executable").getAbsolutePath()
        commandLineDataset = jsap.getString("dataset")  /* default is null, short='d' */
        virtualSchemaName = jsap.getString("virtual-schema-name")   /* default is 'default', short='v' */
        serverPrefix = jsap.getString("server-prefix")        /* default is 'uswest', short is 'u' */
        outputFolder = jsap.getString("output-folder")        /* default is ".", short is 'f' */
        organism = jsap.getString("organism")            /* default is organismToDatasetPrefix.keySet().asList()[0], short is 'r' */
        chromListFile = jsap.getFile("chromosome-list-file")
        if (chromListFile != null) {
            if (!chromListFile.exists()) {
                println "The specified -c/--chromosome-list-file ${chromListFile} didn't exist."
                return false
            }
        }

        if (!organismToDatasetPrefix[organism]) {
            def autoGeneratedDatasetName = buildDefaultDataset(organism)
            organismToDatasetPrefix[organism] = autoGeneratedDatasetName
            System.err.println("Using automatically generated dataset name ${autoGeneratedDatasetName} for unrecognized organism ${organism}.")
            /*System.err.println("Invalid organism specified (${organism})")
            System.err.println("Currently only ${organismToDatasetPrefix.keySet().join(", ")} are supported")
            return false  */
        }

        String exports = jsap.getString("exports")            /* default is exportTypeToDataMap.keySet().asList()[0], short is e */
        if (exports == 'all') {
            exportsList = exportTypeToDataMap.keySet() as String[]
        } else {
            exportsList = exports.split(",")
            for (String export in this.exportsList) {
                if (!exportTypeToDataMap[export]) {
                    System.err.println("Invalid name specified in -e | --exports")
                    System.err.println("Currently only ${exportTypeToDataMap.keySet().join(", ")} are supported")
                    return false
                }
            }
        }


        if (!verifyExecutable("${TABIX_EXEC}")) {
            System.err.println("Please provide the --tabix-executable option. ");
            return false;

        }

        if (!verifyExecutable("${BGZIP_EXEC}")) {
            System.err.println("Please provide the --bgzip-executable option. ");
            return false;
        }
        if (!jsap.success()) {
            jsapConfig.printUsage();

            for (java.util.Iterator errs = jsap.getErrorMessageIterator();
            errs.hasNext();) {
                System.err.println("Error: " + errs.next());
            }

            System.err.println();
            System.err.println("Usage: grooby Biomart.groovy ");


            System.exit(1);
        }
        return jsap.success()

    }

    def buildDefaultDataset(String organism) {
        int lastPartIndex = organism.indexOf("_") + 1
        def result = organism.charAt(0).toString() + organism.subSequence(lastPartIndex, organism.length())
        return result.toString().toLowerCase()
    }

    public boolean verifyExecutable(final String path) {
        File file = new File(path)
        if (file.exists() && file.canRead() && file.canExecute()) {
            return true
        } else {
            println "Executable ${path} isn't available."
            return false
        }
    }
}
