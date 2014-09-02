/*
 * Copyright (c) 2011  by Cornell University and the Cornell Research
 * Foundation, Inc.  All Rights Reserved.
 * 
 * Permission to use, copy, modify and distribute any part of GobyWeb web 
 * application for next-generation sequencing data alignment and analysis, 
 * officially docketed at Cornell as D-5061 ("WORK") and its associated 
 * copyrights for educational, research and non-profit purposes, without 
 * fee, and without a written agreement is hereby granted, provided that 
 * the above copyright notice, this paragraph and the following three 
 * paragraphs appear in all copies.
 * 
 * Those desiring to incorporate WORK into commercial products or use WORK 
 * and its associated copyrights for commercial purposes should contact the 
 * Cornell Center for Technology Enterprise and Commercialization at 
 * 395 Pine Tree Road, Suite 310, Ithaca, NY 14850; 
 * email:cctecconnect@cornell.edu; Tel: 607-254-4698; 
 * FAX: 607-254-5454 for a commercial license.
 * 
 * IN NO EVENT SHALL THE CORNELL RESEARCH FOUNDATION, INC. AND CORNELL 
 * UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF 
 * WORK AND ITS ASSOCIATED COPYRIGHTS, EVEN IF THE CORNELL RESEARCH FOUNDATION, 
 * INC. AND CORNELL UNIVERSITY MAY HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH 
 * DAMAGE.
 * 
 * THE WORK PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE CORNELL RESEARCH 
 * FOUNDATION, INC. AND CORNELL UNIVERSITY HAVE NO OBLIGATION TO PROVIDE 
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE CORNELL 
 * RESEARCH FOUNDATION, INC. AND CORNELL UNIVERSITY MAKE NO REPRESENTATIONS AND 
 * EXTEND NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT
 * NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF WORK AND ITS ASSOCIATED COPYRIGHTS 
 * WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
 */

import com.martiansoftware.jsap.JSAPResult
import edu.cornell.med.icb.util.ICBStringUtils
import edu.cornell.med.icb.goby.modes.ConcatenateCompactReadsMode
import edu.cornell.med.icb.goby.modes.CompactFileStatsMode
import edu.cornell.med.icb.goby.modes.FastaToCompactMode
import edu.cornell.med.icb.goby.modes.SampleQualityScoresMode
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet
import org.apache.commons.io.FileUtils
import org.apache.commons.io.FilenameUtils
import org.apache.commons.logging.Log
import org.apache.commons.logging.LogFactory
import java.util.regex.Pattern
import org.campagnelab.groovySupport.ExecAndRemote
import org.campagnelab.groovySupport.JsapSupport
import org.campagnelab.groovySupport.sample.SampleDetails
import org.campagnelab.groovySupport.sample.SampleDetailsUtil


public class TaskProcessSample {

    private static final Log LOG = LogFactory.getLog(TaskProcessSample.class);

    private final static String JSAP_XML_CONFIG = """
        <jsap>
            <parameters>
                <flaggedOption>
                    <id>output-stats-filename</id>
                    <longFlag>output-stats</longFlag>
                    <required>true</required>
                    <help>The filename were output statistics will be written.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>goby-jar-dir</id>
                    <longFlag>goby-jar-dir</longFlag>
                    <required>true</required>
                    <help>The directory where the goby.jar that will be use exists.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>jvm-flags</id>
                    <longFlag>jvm-flags</longFlag>
                    <required>true</required>
                    <help>Java JVM flags to use.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>cluster-reads-dir</id>
                    <longFlag>cluster-reads-dir</longFlag>
                    <required>true</required>
                    <help>The directory on the cluster where the reads file live.</help>
                </flaggedOption>
                <unflaggedOption>
                    <id>web-sample-files</id>
                    <required>true</required>
                    <greedy>true</greedy>
                    <help>The source sample files (full path, fastq, fastq, compact-reads) on the web server to process for ONE sample. These will be stored on the cluster in --cluster-reads-dir.</help>
                </unflaggedOption>
                <flaggedOption>
                    <id>color-space</id>
                    <longFlag>color-space</longFlag>
                    <required>true</required>
                    <help>Set to true if the sample is colorspace.</help>
                    <stringParser>
                        <classname>BooleanStringParser</classname>
                    </stringParser>
                </flaggedOption>
                <flaggedOption>
                    <id>sample-tag</id>
                    <longFlag>sample-tag</longFlag>
                    <required>true</required>
                    <help>The tag for the sample for which we are processing reads.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>sample-name</id>
                    <longFlag>sample-name</longFlag>
                    <required>true</required>
                    <help>The sample name, should be filename friendly (no spaces, etc.).</help>
                </flaggedOption>
                <flaggedOption>
                    <id>merge-plan-filename</id>
                    <longFlag>merge-plan-filename</longFlag>
                    <required>true</required>
                    <help>The name of the file that contains the merge plan (formerly details.tsv file).</help>
                </flaggedOption>
                <flaggedOption>
                    <id>platform</id>
                    <longFlag>platform</longFlag>
                    <required>true</required>
                    <help>The platform for the sample.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>quality-encoding</id>
                    <longFlag>quality-encoding</longFlag>
                    <required>true</required>
                    <help>The quality encoding for the reads, in the case the reads are fastq.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>first-file-tag</id>
                    <longFlag>first-file-tag</longFlag>
                    <required>true</required>
                    <help>The tag of the first file in the list of files to be processed, if there are more than one we are performing a merge.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>web-files-dir</id>
                    <longFlag>web-files-dir</longFlag>
                    <required>true</required>
                    <help>The destination directory on the web server machine for returning stats files, etc.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>ssh-prefix</id>
                    <longFlag>ssh-prefix</longFlag>
                    <required>true</required>
                    <help>The ssh / scp prefix (such as username@server).</help>
                </flaggedOption>
                <flaggedOption>
                    <id>queue-writer-prefix</id>
                    <longFlag>queue-writer-prefix</longFlag>
                    <required>false</required>
                    <help>The command line prefix for writing to the queue</help>
                </flaggedOption>
                <flaggedOption>
                     <id>queue-writer-prefix-variable</id>
                     <longFlag>queue-writer-prefix-variable</longFlag>
                     <required>false</required>
                     <help>The name of an environment variable containing the command line prefix for writing to the queue.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>job-start-status</id>
                    <longFlag>job-start-status</longFlag>
                    <required>false</required>
                    <help>The generic "start" status message to use when writing to the queue.</help>
                </flaggedOption>
                <flaggedOption>
                    <id>work-dir</id>
                    <longFlag>work-dir</longFlag>
                    <required>false</required>
                    <stringParser>
                        <classname>StringStringParser</classname>
                    </stringParser>
                    <defaults>
                        <string>/tmp</string>
                    </defaults>
                    <help>A work / temporary directory to work it. It is assumed that any files
                          this process creates in this directory are no longer needed after
                          this process quotes and will be cleaned up by a separate process
                          (such as Oracle Grid Engine).</help>
                </flaggedOption>
                <switch>
                    <id>skip-queue-writing</id>
                    <longFlag>skip-queue-writing</longFlag>
                    <help>If this flag is enabled, the GobyWeb queue will not be written to.</help>
                </switch>
                <switch>
                    <id>help</id>
                    <longFlag>help</longFlag>
                    <help>Help.</help>
                </switch>
            </parameters>
        </jsap>
    """

    public final static COMPACT_READS_ALLOWED_EXTS = [".compact-reads"]
    public final static FASTX_ALLOWED_EXTS = [
            ".fa", ".fq", ".fasta", ".fastq", ".csfasta",
            ".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz", ".csfasta.gz", ".txt.gz",
            ".fastq.gz.tar"]
    public final static ALL_ALLOWED_EXTS = []
    static {
        ALL_ALLOWED_EXTS.addAll COMPACT_READS_ALLOWED_EXTS
        ALL_ALLOWED_EXTS.addAll FASTX_ALLOWED_EXTS
    }

    /**
     * For fastq.gz.tar files, this will pull the basename of the file
     * plus the pair indicator ("_R1_", "_R2_" in the filename will return "1", "2").
     */
    private final Pattern TAR_FILENAME_PAIRS_PATTERN = ~/^(.*)_R([12]?)_[0-9]+\.fastq\.gz$/

    // Configuration that would be fed in from cluster script
    String gobyJarDir
    String jvmFlags
    String clusterReadsDir
    String[] webSampleFiles
    String sampleTag
    String firstFileTag
    String qualityEncoding
    String platform
    String sampleName
    boolean colorSpace
    String sshPrefix
    String webFilesDir
    String webUploadsDir
    String queueWriterPrefix
    String queueJobStartStatusCode = "UNKNOWN"
    String workDir = "/tmp"
    String outStatsFilename
    // This will be populated during execution
    int numberOfReads
    Writer statsWriter
    String readQualChartFilename
    String mergePlanFilename
    ExecAndRemote exec = new ExecAndRemote()

    public static void main(String[] args) {
        def exec = new TaskProcessSample()
        if (args.length == 1 && args[0] == "--test") {
            exec.testConfigure()
        } else {
            exec.configure(args)
        }
        int retval = exec.processSingleSample()
        if (retval != 0) {
            System.exit retval
        }
    }

    public TaskProcessSample() {
        // This will be populated during execution
        numberOfReads = -1
        statsWriter = null
        readQualChartFilename = null
        exec = new ExecAndRemote()
    }

    private void testConfigure() {
        gobyJarDir =null
        jvmFlags = "-Xmx2g"
        clusterReadsDir = "/home/gobyweb/GOBYWEB_FILES-kdorff/kdorff"
        webSampleFiles = [
                "/home/gobyweb/GOBYWEB_UPLOADS/01.compact-reads",
                "/home/gobyweb/GOBYWEB_UPLOADS/02.compact-reads",
                "/home/gobyweb/GOBYWEB_UPLOADS/03.compact-reads",
                "/home/gobyweb/GOBYWEB_UPLOADS/04.compact-reads"] as String[]
        sampleTag = "KEWEFWD"
        firstFileTag = "JGEWKQS"
        qualityEncoding = "Illumina"
        platform = "Illumina"
        sampleName = "kevinreads"
        colorSpace = false
        sshPrefix = "kdorff@mac133990.med.cornell.edu"
        webFilesDir = "/Users/gobyweb/GOBYWEB_FILES-kdorffmac/kdorff"
        exec.setSkipQueueWriting(true)
    }

    private void configure(String[] args) {
        final JSAPResult jsapResult = new JsapSupport()
                .setScriptName("TaskProcessSample.groovy")
                .setArgs(args)
                .setXmlConfig(JSAP_XML_CONFIG)
                .parse()
        webSampleFiles = jsapResult.getStringArray("web-sample-files")
        gobyJarDir = jsapResult.getString("goby-jar-dir")
        jvmFlags = jsapResult.getString("jvm-flags")
        mergePlanFilename = jsapResult.getString("merge-plan-filename")
        clusterReadsDir = jsapResult.getString("cluster-reads-dir")
        sampleTag = jsapResult.getString("sample-tag")
        firstFileTag = jsapResult.getString("first-file-tag")
        qualityEncoding = jsapResult.getString("quality-encoding")
        platform = jsapResult.getString("platform")
        sampleName = jsapResult.getString("sample-name")
        colorSpace = jsapResult.getBoolean("color-space")
        sshPrefix = jsapResult.getString("ssh-prefix")
        webFilesDir = jsapResult.getString("web-files-dir")
        queueWriterPrefix = jsapResult.getString("queue-writer-prefix")
        outStatsFilename = jsapResult.getString("output-stats-filename")

        String queueWriterPrefixVariable = jsapResult.getString("queue-writer-prefix-variable")

        if (queueWriterPrefixVariable != null && queueWriterPrefix == null) {
            queueWriterPrefix = System.getenv(queueWriterPrefixVariable)
        }
        queueJobStartStatusCode = jsapResult.getString("job-start-status")
        workDir = jsapResult.getString("work-dir")
        if (queueWriterPrefix && queueJobStartStatusCode) {
            exec.setupQueue(queueWriterPrefix, queueJobStartStatusCode)
            exec.setSkipQueueWriting(jsapResult.getBoolean("skip-queue-writing"))
        } else {
            exec.setSkipQueueWriting(true)
        }
    }

    /**
     * Process the uploaded samples.

     */
    private int processSingleSample() {
        int retval = 99
        println "Processing sample.tag=${sampleTag}"

        // Copy files from web server
        List<String> processFilePaths = new ArrayList<String>()
        int i = 0
        int numFiles = webSampleFiles.length
        final String destinationDir = "CONVERTED"
        for (String webSampleFile in webSampleFiles) {

            String localFilename = FilenameUtils.getName(webSampleFile)

            // copy all files to CONVERTED
            // if input file is a compact-reads, we copy it to the CONVERTED folder:
            if (FilenameUtils.getExtension(webSampleFile).equals("compact-reads")) {
                def destinationFile = new File(FilenameUtils.concat(destinationDir, localFilename))

                println "Copying ${webSampleFile} to ${destinationFile}"

                FileUtils.copyFile(new File(webSampleFile), destinationFile)
                processFilePaths.add(destinationFile.getAbsolutePath())
            } else {
                processFilePaths.add(new File(webSampleFile).getAbsolutePath())
            }
        }

        readQualChartFilename = "${clusterReadsDir}/${firstFileTag}.quality-stats.tsv"

        final boolean pairedSamplesFound
        final String sampleName = null
        final Map<String, SampleDetails> filenameToSampleDetailsMap

        // When the detailsFile is read from here, it should contain files for ONE sample name
        final Map<String, List<SampleDetails>> detailsMap = SampleDetailsUtil.readSampleDetails(mergePlanFilename)
        println "Map before removing pairs: ${detailsMap}"
        pairedSamplesFound = SampleDetailsUtil.removePairs(detailsMap)
        sampleName = (detailsMap.keySet() as List)[0]
        filenameToSampleDetailsMap = SampleDetailsUtil.sampleDetailsMapToFilenameLcMap(detailsMap)
        if (pairedSamplesFound) {
            println "Map after removing pairs: ${detailsMap}"
            // Since we've removed pairs, we can re-obtain processFilenames from the detailsMap.
            // Make a map of the process filenames without the tag to the process filename
            Map<String, String> noTagToProcessMap = [:]
            processFilePaths.each { String filePath ->
                noTagToProcessMap[new File(filePath).getName().substring(8)] = filePath
            }
            // Copy just the ones that aren't part 1 of the pair (just part 0)
            processFilePaths.clear()
            (detailsMap[sampleName])*.filename.each { final String filename ->
                processFilePaths.add(noTagToProcessMap[filename])
            }
            renamePairFiles(filenameToSampleDetailsMap, processFilePaths, noTagToProcessMap)
        } else {
            println "No pairs found in map"
        }
        println "filenameToSampleDetailsMap=${filenameToSampleDetailsMap}"

        println "processFilePaths=${processFilePaths}"


        // Convert any fasta/fastq to compact-reads
        final List<File> stepTwoFiles = new ArrayList<File>()

        String pairFilename
        for (String processFilePath in processFilePaths) {
            pairFilename = null
            if (processFilePath.endsWith(".compact-reads")) {

                stepTwoFiles << new File(processFilePath)
            } else {
                boolean processFqGzTar = false
                if (processFilePath.endsWith(".fastq.gz.tar")) {
                    processFqGzTar = true
                    (processFilePath, pairFilename) = processFastqGzTar(workDir, processFilePath)
                    if (processFilePath == null) {
                        retval = 4
                        return retval
                    }
                }
                exec.queueMessage sampleTag, "Converting reads from fasta/fastq to Goby Compact-reads format"
                String fafqFilename = "${processFilePath}"
                String localFafqFilename = processFilePath

                String outputBasename
                if (processFqGzTar) {
                    outputBasename = removeFileExtension(fafqFilename)
                    if (pairFilename && outputBasename.endsWith("_1")) {
                        // Remove the pair designation on the output compact-reads file
                        outputBasename = outputBasename[0..-3]
                    }
                } else {
                    if (processFilePath.size() > 1) {
                        outputBasename = removeFileExtension(fafqFilename)
                    } else {
                        outputBasename = sampleName ?: removeFileExtension(fafqFilename)
                    }
                }

                String compactFilename = "${outputBasename}.compact-reads"
                def destinationFile = new File(FilenameUtils.concat(destinationDir, compactFilename))

                String localCompactFilename = destinationFile.getAbsolutePath()
                try {

                    // This will try to automatically set the quality encoding based
                    // on the results of SampleQualityScoresMode.
                    SampleQualityScoresMode sqs = new SampleQualityScoresMode()
                    sqs.addInputFilename localFafqFilename
                    sqs.execute()
                    List<String> likelyEncodingList = sqs.likelyEncodings
                    String likelyEncoding
                    if (likelyEncodingList) {
                        likelyEncoding = likelyEncodingList[0]
                        if (likelyEncoding == "Illumina/Solexa") {
                            likelyEncoding = "Illumina"
                        }
                        likelyEncoding = likelyEncoding.toUpperCase()
                    } else {
                        println "Determination of Likely encoding failed."
                    }

                    FastaToCompactMode convert = new FastaToCompactMode()
                    convert.addInputFilename localFafqFilename
                    convert.setOutputFilename localCompactFilename
                    if (processFqGzTar) {
                        if (pairFilename) {
                            convert.processPairs = true
                            convert.pairIndicator1 = "_1"
                            convert.pairIndicator2 = "_2"
                        }
                    } else {
                        if (pairedSamplesFound) {
                            // Original filename lowercase without tag from processFilename
                            final String originalFilenameLc = processFilePath.toLowerCase().substring(8)
                            final SampleDetails sampleDetails = filenameToSampleDetailsMap[originalFilenameLc]
                            def pairIndicatorsList = sampleDetails?.pairIndicatorsList()
                            if (pairIndicatorsList) {
                                convert.processPairs = true
                                convert.pairIndicator1 = pairIndicatorsList[0]
                                convert.pairIndicator2 = pairIndicatorsList[1]
                            } else {
                                println "WARNING: It was believed the samples were paired but no pair indicators were found."
                            }
                        }
                    }
                    println "Convert ${localFafqFilename} to ${localCompactFilename} pairedSamplesFound=${pairedSamplesFound} likelyEncoding=${likelyEncoding}"
                    if (pairedSamplesFound) {
                        println ".. convert.processPairs=${convert.processPairs}"
                        println ".. convert.pairIndicator1=${convert.pairIndicator1}"
                        println ".. convert.pairIndicator2=${convert.pairIndicator2}"
                    }

                    // Changed from manual selection of encoding to automatic selection.
                    if (likelyEncoding) {
                        if (likelyEncoding != "FASTA") {
                            convert.qualityEncoding = likelyEncoding
                        }
                    } else {
                        if (qualityEncoding != "Other") {
                            convert.setQualityEncoding qualityEncoding
                        }
                    }

                    println "Conversion starting"
                    convert.execute()
                    println "Conversion finished"
                    stepTwoFiles << new File(localCompactFilename)

                } catch (IOException e) {
                    // The conversion failed. Nothing else can be done.
                    println "Conversion to compact-reads failed, IOException ${e.message}"
                    exec.queueMessage sampleTag, "Conversion to compact-reads failed"
                    return
                } catch (IllegalArgumentException e) {
                    // The conversion failed. Nothing else can be done.
                    println "Conversion to compact-reads failed, IllegalArgumentException ${e.message}"
                    exec.queueMessage sampleTag, "Conversion to compact-reads failed"
                    return
                }
            }
        }

        try {
            // TODO: if we are concat'ing, we should probaly check the first read-index of
            // TODO: each of the input files to make sure we don't duplicate.
            println "stepTwoFiles=${stepTwoFiles}"
            File localFile
            final String storedName =
                    "${firstFileTag}-${ICBStringUtils.safeFilename(FilenameUtils.getBaseName(sampleName))}.compact-reads"
            if (stepTwoFiles.size() > 1) {
                exec.queueMessage sampleTag, "Concatenating reads"
                ConcatenateCompactReadsMode concat = new ConcatenateCompactReadsMode();
                localFile = new File(storedName)
                concat.setOutputFilename localFile.toString()
                concat.setQuickConcat false
                for (stepTwoFile in stepTwoFiles) {
                    println "Adding file to concat ${stepTwoFile}"
                    concat.addInputFile stepTwoFile
                }
                concat.execute()

                for (stepTwoFile in stepTwoFiles) {
                    FileUtils.deleteQuietly(stepTwoFile)
                }
            } else {
                File prevLocalFile = stepTwoFiles[0]
                localFile = new File("${storedName}")
                localFile.delete() //delete target if exists
                println "Single file. renaming ${prevLocalFile} to ${localFile}"
                FileUtils.copyFile(prevLocalFile, localFile)
            }
            println "Preparing stats file..."
            //String mergePlanFilename = "${clusterReadsDir}/${firstFileTag}.details.txt"
            statsWriter = (new File(outStatsFilename)).newWriter(false)
            statsWriter.writeLine "ngFile.tag=${firstFileTag}"
            statsWriter.writeLine "ngFile.storedName=${storedName}"
            statsWriter.writeLine "ngFile.storedDir=${clusterReadsDir}"
            statsWriter.writeLine "ngFile.size=${localFile.length()}"
            if (localFile && localFile.exists()) {
                println "Generating reads length..."
                getReadLengths(localFile)
                getReadQualityStats(localFile)
                calculateHeptamersWeights(localFile)
                if (platform != "SOLiD") {
                    calculateGCWeights(localFile)
                }
            } else {
                println "local file does not exist..."
                retval = 2
            }
            statsWriter.writeLine "sample.readyToAlign=true"
            statsWriter.close()
            statsWriter = null
            //}
            retval = 0
        } catch (Exception e) {
            println "Error doing background processing of Sample with sample.tag=$sampleTag ${e.message}"
            e.printStackTrace()
            exec.queueMessage sampleTag, "Error processing sample"
            retval = 1
        } finally {
            // If failed, notify via queue
            if (statsWriter != null) {
                statsWriter.writeLine "sample.readyToAlign=false"
                statsWriter.close()
            }

        }
        return retval
    }

    /**
     * Return two lists of filenames. The first list is the R1 files (first set of reads
     * if paired end). The second list will be empty unless the reads are paired in,
     * in which case the list of filenames that are the pairs for the first list.
     * The order of the files in the two lists should be matched.
     */
    def tarFilenames(dir, srcFilename) {
        def (exit, filenames) = exec.execReturnStdoutList("tar -tf ${dir}/${srcFilename}")
        // Sort so the the filenames so if we need to pair them up, they will be in the right order
        Collections.sort(filenames)
        // These files will be exactracted to workDir, so prepend this to each filename.
        filenames = filenames.collect { String filename ->
            "${workDir}/${filename}".toString()
        }

        def pairFilenames = []
        def removeFilenames = []
        for (final String filename in filenames) {
            if (filename.endsWith("/")) {
                removeFilenames << filename
            } else {
                def noPathFilename = FilenameUtils.getName(filename)
                noPathFilename.find(TAR_FILENAME_PAIRS_PATTERN) { match, basename, pairIndicator ->
                    if (pairIndicator == "2") {
                        pairFilenames << filename
                    }
                }
            }
        }
        for (final String removeFilename in removeFilenames) {
            filenames.remove(removeFilename)
        }
        for (final String pairFilename in pairFilenames) {
            filenames.remove(pairFilename)
        }
        return [filenames, pairFilenames]
    }

    /**
     * Input file is something.fastq.gz.tar
     * Create something.fastq.gz from it.
     * This method will delete the fastq.gz.tar file after it has been processed
     * @param srcFilename the full path to the source filename
     * @return new filename something.fastq.gz
     */
    private List<String> processFastqGzTar(final String dir, final String srcFilename) {
        if (!srcFilename.endsWith(".tar")) {
            return [srcFilename, null]
        }
        exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz"

        List<String> filenames, pairFilenames
        (filenames, pairFilenames) = tarFilenames(dir, srcFilename)
        if (!filenames) {
            println "No filenames found in tar file"
            return [null, null]
        }
        if (pairFilenames && filenames.size() != pairFilenames.size()) {
            println "Pair files found, but the number of filenames and pairFilenames don't match."
        }
        boolean pairedEnd = pairFilenames.size() > 0

        String destFilename
        String pairFilename
        if (!pairedEnd) {
            // Strip off .tar
            destFilename = srcFilename[0..-5]
        } else {
            // Strip off .fastq.gz.tar then add _1/_2 suffix and .fastq.gz back on
            destFilename = srcFilename[0..-14] + "_1.fastq.gz"
            pairFilename = srcFilename[0..-14] + "_2.fastq.gz"
        }

        // Non-paired end. This is the fast method as the files don't need to be extracted
        // to create the fastq file.
        if (!pairedEnd) {
            int retval = exec.bashExec("tar -xf ${dir}/${srcFilename} --wildcards --no-anchored '*_R1_*' -O | gunzip | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--\$' | gzip > ${dir}/${destFilename}")
            if (retval == 0) {
                exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz completed"
            } else {
                exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz failed"
                return [null, null]
            }
        } else {
            // As tar files are not random accesss
            // and the order of the files within the tar are completely random,
            // we need to extract all of the tar files to get to the files to filter/merge
            // together in the correct order.

            // Extract the files from the tar
            int retval = exec.bashExec("tar -xf ${dir}/${srcFilename} -C ${workDir}")
            if (retval != 0) {
                exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz failed, error extracting files form tar"
                return [null, null]
            }
            retval = exec.bashExec("zcat ${filenames.join(" ")} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--\$' | gzip > ${dir}/${destFilename}")
            if (retval != 0) {
                exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz failed, error concat/filter R1 reads."
                return [null, null]
            }
            retval = exec.bashExec("zcat ${pairFilenames.join(" ")} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--\$' | gzip > ${dir}/${pairFilename}")
            if (retval != 0) {
                exec.queueMessage sampleTag, "Processing .fastq.gz.tar to .fastq.gz failed, error concat/filter R2 reads."
                return [null, null]
            }
        }

        return [destFilename, pairFilename]
    }

    /**
     * Assuming the extension is listed in ALL_ALLOWED_EXTS this will remove
     * the file extension.
     */
    private String removeFileExtension(String filename) {
        if (!filename) {
            return filename
        }
        String foundExt = null
        for (goodExt in ALL_ALLOWED_EXTS) {
            if (filename.endsWith(goodExt)) {
                foundExt = goodExt
                break
            }
        }
        if (!foundExt) {
            return filename
        }
        return filename.substring(0, filename.size() - foundExt.size())
    }

    private boolean getReadLengths(final File sample) {
        exec.queueMessage sampleTag, "Determining read lengths"
        println "getReadLengths for ${sample.getAbsolutePath()}"
        CompactFileStatsMode stats = new CompactFileStatsMode();
        if (sample && sample.exists()) {
            println "Adding file to stats ${sample}"
            stats.addInputFile sample
            println "Fetching stats for file(s)"
            stats.execute()
            numberOfReads = stats.numberOfReads
            statsWriter.writeLine "sample.minimumLength=${stats.minReadLength}"
            statsWriter.writeLine "sample.maximumLength=${stats.maxReadLength}"
            statsWriter.writeLine "sample.numberOfReads=${stats.numberOfReads}"
            statsWriter.writeLine "sample.pairedSample=${stats.allPairedSamples}"
            return true
        } else {
            println "getReadLengths: sample does not exist"
            return false
        }
    }

    private void calculateHeptamersWeights(final File sample) {
        def completed = false
        try {
            def heptamerInfoFile = new File("${gobyJarDir}/heptamer-${platform}-info.bin")
            if (sample && sample.exists() && heptamerInfoFile.exists()) {
                String colorSpaceOpt = ""
                if (colorSpace) {
                    colorSpaceOpt = "--color-space"
                }
                String cl = "java ${jvmFlags} -jar ${gobyJarDir}/goby.jar -m reads-to-weights ${sample} ${colorSpaceOpt} --method heptamers --heptamer-info ${heptamerInfoFile}"
                println "Executing read-quality-stats / heptamers mode for ${sample}"
                exec.queueMessage sampleTag, "Creating heptamer weights file"
                def exitValue = exec.exec(cl)
                if (exitValue == 0) {
                    println "... done!"
                    completed = true
                } else {
                    println "Failed with exitValue=${exitValue}"
                }
            }
        } catch (Exception e) {
            LOG.error("Error creating heptamers weights", e)
        } finally {
            statsWriter.writeLine "sample.attributes.heptamersWeightsAvailable=${completed}"
        }
    }

    private void calculateGCWeights(final File sample) {
        def completed = false
        try {
            if (sample && sample.exists()) {
                String colorSpaceOpt = ""
                if (colorSpace) {
                    colorSpaceOpt = "--color-space"
                }
                String cl = "java ${jvmFlags} -jar ${gobyJarDir}/goby.jar -m reads-to-weights ${sample} ${colorSpaceOpt} --method gc "
                println "Executing read-quality-stats / gc mode for ${sample}"
                exec.queueMessage sampleTag, "Creating GC weights file"
                def exitValue = exec.exec(cl)
                if (exitValue == 0) {
                    println "... done!"
                    completed = true
                } else {
                    println "Failed with exitValue=${exitValue}"
                }
            }
        } catch (Exception e) {
            LOG.error("Error creating GC weights", e)
        } finally {
            statsWriter.writeLine "sample.attributes.gcWeightsAvailable=${completed}"
        }
    }

    private void getReadQualityStats(final File sample) {
        // TODO: add to queue "Determining read quality"
        println "getReadQualityStats for ${sample.getAbsolutePath()}"
        if (sample && sample.exists()) {
            String cl = "java ${jvmFlags} -jar ${gobyJarDir}/goby.jar -m read-quality-stats -o ${readQualChartFilename} ${sample} "
            if (numberOfReads < 5000) {
                cl += "-p 1.0 "
            }
            println "Executing read-quality-stats mode"
            exec.queueMessage sampleTag, "Determining reads quality"
            def exitValue = exec.exec(cl)
            if (exitValue == 0) {
                println "... done!"
            } else {
                println "Failed with exitValue=${exitValue}"
            }
        } else {
            println "Could not find sample file ${sample} to process"
        }
    }

    /**
     * Goby's fasta-to-compact expects the two pair filenames to be exactly the same except for pair
     * indicators. In actuality here, the tags will be different. This will rename the pair file (with
     * the second file in the pair) to have the same tag as the first file in the pair
     * @param filenameToSampleDetailsMap
     * @param processFilenames
     * @param noTagToProcessMap
     */
    private renamePairFiles(final Map<String, SampleDetails> filenameToSampleDetailsMap,
                            final List<String> processFilenames,
                            final Map<String, String> noTagToProcessMap) {
        for (final String processFilename in processFilenames) {
            final String firstInPairNoTagLc = processFilename.toLowerCase().substring(8).toLowerCase()
            final SampleDetails sampleDetails = filenameToSampleDetailsMap[firstInPairNoTagLc]
            if (sampleDetails) {
                final String pairKey = makePairFilename(firstInPairNoTagLc, sampleDetails.pairIndicatorsList())
                if (pairKey) {
                    def originalPairFilename = noTagToProcessMap[pairKey]
                    if (originalPairFilename) {
                        final String newPairFilename = makePairFilename(processFilename, sampleDetails.pairIndicatorsList())
                        final File renameFrom = new File("${workDir}/${originalPairFilename}")
                        final File renameTo = new File("${workDir}/${newPairFilename}")
                        println "Renaming pair from ${renameFrom} to ${renameTo} (first is ${processFilename})"
                        renameFrom.renameTo(renameTo)
                    }
                } else {
                    println "pairKey=${pairKey} not found"
                }
            } else {
                println "sampleDetails for ${firstInPairNoTagLc} not found"
            }
        }
    }

    /**
     * Give the filename and the pair indicators, provide a filename which switches the pair
     * indicator. It is assumed that the pair indicator is the LAST occurence of that string
     * in the filename).
     * @param firstFilename
     * @param pairIndicators
     * @return
     */
    public makePairFilename(final String firstFilename, final List<String> pairIndicators) {
        if (!firstFilename || !pairIndicators) {
            return null
        }
        final String pairIndicator1 = pairIndicators[0]
        final String pairIndicator2 = pairIndicators[1]
        final int firstTokenPos = firstFilename.lastIndexOf(pairIndicator1);
        if (firstTokenPos == -1) {
            // No pairIndicator1 token. Not the first file in a pair.
            return null;
        }
        final StringBuilder pairFilenameSb = new StringBuilder();
        pairFilenameSb.append(firstFilename.substring(0, firstTokenPos));
        pairFilenameSb.append(pairIndicator2);
        pairFilenameSb.append(firstFilename.substring(firstTokenPos + pairIndicator1.length()));
        return pairFilenameSb.toString();
    }
}

