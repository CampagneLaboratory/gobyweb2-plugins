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
import edu.cornell.med.icb.goby.readers.vcf.ColumnType
import edu.cornell.med.icb.goby.readers.vcf.VCFParser
import groovy.sql.Sql
import java.sql.SQLException
import org.apache.commons.io.FilenameUtils;
import it.unimi.dsi.logging.ProgressLogger
import org.apache.log4j.Level
import org.apache.log4j.Logger
import java.sql.PreparedStatement
import edu.cornell.med.icb.goby.readers.vcf.GroupAssociations
import org.campagnelab.groovySupport.Sqlite
import org.campagnelab.groovySupport.JsapSupport
import org.campagnelab.groovySupport.ExecAndRemote
import org.campagnelab.lucene.analysis.SpecialSymbolAnalyzer
import org.apache.lucene.document.NumericField
import org.apache.lucene.document.Field
import org.apache.lucene.document.Document
import org.apache.lucene.index.IndexWriter
import org.apache.lucene.index.IndexWriterConfig
import org.apache.lucene.store.SimpleFSDirectory
import org.apache.lucene.util.Version
import org.apache.commons.io.FileUtils
import org.campagnelab.groovySupport.db.DbColumnDetails
import org.apache.lucene.document.Fieldable

/**
 * Create a Sqlite database from a TSV or VCF file that came from a GobyWeb Alignment Comparison (DiffExp).
 * This will be run on the server as the last step of any execution that creates vcf/tsv output that needs
 * to be browsed.
 */
public class TsvVcfToSqlite {

    private static final Logger LOG = Logger.getLogger(TsvVcfToSqlite.class);

    final boolean DEBUG = false

    private static final String JSAP_XML_CONFIG = """
        <jsap>
            <parameters>
                <flaggedOption>
                    <id>tag</id>
                    <longFlag>tag</longFlag>
                    <shortFlag>t</shortFlag>
                    <required>false</required>
                    <help>The GobyWeb tag for the file to process. If not provided, it will assume the start of the filename before the first "-" is the tag. If more than one input file is provided, this will be ignored.</help>
                </flaggedOption>
                <unflaggedOption>
                    <id>input</id>
                    <required>true</required>
                    <greedy>true</greedy>
                    <help>The tsv/vcf file(s) to convert to Sqlite.</help>
                </unflaggedOption>
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
                    <help>The name of an environment variable containing the command line prefix for writing to the queue</help>
                </flaggedOption>
                <flaggedOption>
                    <id>job-start-status</id>
                    <longFlag>job-start-status</longFlag>
                    <required>false</required>
                    <help>The generic "start" status message to use when writing to the queue</help>
                </flaggedOption>
                <flaggedOption>
                    <id>export-format</id>
                    <longFlag>export-format</longFlag>
                    <shortFlag>f</shortFlag>
                    <required>false</required>
                    <help>The export format, 'sqlite' or 'lucene'.</help>
                    <defaults>
                        <string>sqlite</string>
                    </defaults>
                </flaggedOption>
                <flaggedOption>
                    <id>number-to-export</id>
                    <longFlag>number-to-export</longFlag>
                    <shortFlag>n</shortFlag>
                    <required>false</required>
                    <stringParser>
                        <classname>LongStringParser</classname>
                    </stringParser>
                    <defaults>
                        <string>0</string>
                    </defaults>
                    <help>The requested maximum number of rows to export. A value less than or equal to 0 indicates to export all rows.</help>
                </flaggedOption>
                <switch>
                    <id>console-logging</id>
                    <longFlag>console-logging</longFlag>
                    <shortFlag>l</shortFlag>
                    <help>Enable console logging output without the need to otherwise configure logging.</help>
                </switch>
                <switch>
                    <id>help</id>
                    <longFlag>help</longFlag>
                    <help>Help</help>
                </switch>
                <switch>
                    <id>skip-queue-writing</id>
                    <longFlag>skip-queue-writing</longFlag>
                    <help>If this flag is enabled, the GobyWeb queue will not be written to</help>
                </switch>
            </parameters>
        </jsap>
    """

    // Configured from command line
    String tag
    Set<String> inputFilenames
    String queueWriterPrefix
    String queueJobStartStatusCode = "UNKNOWN"
    long numberRowsToExport = 0

    // Non command line configuration
    ExecAndRemote exec = new ExecAndRemote()

    TsvVcfExporters converter = TsvVcfExporters.SqliteExporter

    public static void main(String[] args) {
        def tsvImporter = new TsvVcfToSqlite()
        if (args.length == 1 && args[0] == "--test") {
            tsvImporter.testConfigure()
        } else {
            tsvImporter.configure(args)
        }
        int retval = tsvImporter.execute()
        if (retval != 0) {
            System.exit retval
        }
    }

    public static void setupConsoleLogging() {
        Logger rootLogger = Logger.getRootLogger();
        rootLogger.setLevel(Level.INFO);
        //rootLogger.addAppender(new ConsoleAppender(new PatternLayout("%-5p [%t]: %m%n")));
    }

    private void testConfigure() {
        // Worked.
        // tag = "RDOTEWM"
        // addInputFilename("/tmp/RDOTEWM-zymo-human-test.sequence-variation-stats.tsv")

        // BAD?
        // tag = "NAJPBOJ"
        // addInputFilename("/tmp/NAJPBOJ.vcf.gz")

        tag = "YVPGOTJ"
        addInputFilename("/tmp/YVPGOTJ.vcf.gz")
        exec.setSkipQueueWriting(true)
    }

    private void configure(String[] args) {
        final JSAPResult jsapResult = new JsapSupport()
                .setScriptName("TsvVcfToSqlite.groovy")
                .setArgs(args)
                .setXmlConfig(JSAP_XML_CONFIG)
                .parse()

        tag = jsapResult.getString("tag")
        addInputFilenames(jsapResult.getStringArray("input"))
        if (inputFilenames.size() > 1) {
            tag = null
        }
        queueWriterPrefix = jsapResult.getString("queue-writer-prefix")
        String queueWriterPrefixVariable = jsapResult.getString("queue-writer-prefix-variable")

        if (queueWriterPrefixVariable != null && queueWriterPrefix == null) {
            queueWriterPrefix = System.getenv(queueWriterPrefixVariable)
        }

        queueJobStartStatusCode = jsapResult.getString("job-start-status")
        numberRowsToExport = jsapResult.getLong("number-to-export")
        if (numberRowsToExport < 0) {
            numberRowsToExport = 0
        }
        if (queueWriterPrefix && queueJobStartStatusCode) {
            exec.setupQueue(queueWriterPrefix, queueJobStartStatusCode)
            exec.setSkipQueueWriting(jsapResult.getBoolean("skip-queue-writing"))
        } else {
            exec.setSkipQueueWriting(true)
        }
        final String converterStr = jsapResult.getString("export-format")
        if (converterStr == "lucene") {
            converter = TsvVcfExporters.LuceneExporter
        } else {
            converter = TsvVcfExporters.SqliteExporter
        }
        if (jsapResult.getBoolean("console-logging")) {
            setupConsoleLogging()
        }
    }

    TsvVcfExporters getConverter() {
        return converter
    }

    void setConverter(final TsvVcfExporters converter) {
        this.converter = converter
    }

    public void addInputFilename(final String inputFilename) {
        if (inputFilenames == null) {
            inputFilenames = new LinkedHashSet<String>()
        }
        if (new File(inputFilename).exists()) {
            inputFilenames << inputFilename
        }
    }

    public void addInputFilenames(final String[] inputFilenames) {
        for (final String inputFilename in inputFilenames) {
            addInputFilename(inputFilename)
        }
    }

    private int execute() {
        String currentTag
        final TsvVcfExporter exporter


        for (final String inputFilename : inputFilenames) {
            // create a new exporter for each file:
            exporter = resetConverter()
            if (tag == null) {
                String subFilename = FilenameUtils.getName(inputFilename)
                String[] parts = subFilename.split("[-\\.]")
                currentTag = parts[0]
            } else {
                currentTag = tag
            }

            final File inputFile = new File(inputFilename)
            if (!inputFile.exists()) {
                println "TSV File ${tsvFile} didn't exist. Not importing."
                LOG.error "TSV File ${tsvFile} didn't exist. Not importing."
                return 1
            }

            final VCFParser parser = new VCFParser(inputFile.getCanonicalPath());
            parser.readHeader();

            final int retval = exporter.convert(inputFile, parser, currentTag)
            if (retval != 0) {
                return retval
            }
        }
        return 0
    }

    private TsvVcfExporter resetConverter() {
        TsvVcfExporter exporter
        if (converter == TsvVcfExporters.SqliteExporter) {
            exporter = new SqliteExporter()
        } else {
            exporter = new LuceneExporter()
        }
        exporter.initialize(exec, numberRowsToExport)
        exporter.DEBUG = DEBUG
        exporter
    }
}

enum TsvVcfExporters {
    SqliteExporter,
    LuceneExporter
}

abstract class TsvVcfExporter {

    boolean DEBUG = false

    ExecAndRemote exec

    File backupFile
    File outputFile

    // the last output filename we imported to

    String outputFilename

    /** Maximum number of rows to export, 0 means all. */
    long numberRowsToExport

    String getOutputFilename() {
        return outputFilename
    }

    void initialize(final ExecAndRemote exec, final long numberRowsToExport) {
        this.exec = exec
        this.numberRowsToExport = numberRowsToExport
    }

    abstract int convert(final File inputFile, final VCFParser parser, final String currentTag)

    def dumpList(list) {
        int i = 0
        println "["
        list.each { item ->
            println "     ${i}:${item}"
            i++
        }
        println "]"
    }

    public void backupExisting() {
        if (outputFile.exists()) {
            def now = new Date()
            String dateStr = "${String.format("%tY%<tm%<td", now)}-${String.format("%tH%<tM%<tS", now)}"
            backupFile = new File("${outputFile.toString()}-${dateStr}.backup")
            println "Backing up ${outputFile} to ${backupFile}"
            outputFile.renameTo(backupFile)
        } else {
            backupFile = null
        }
    }

    public void restoreBackup() {
        if (backupFile != null && backupFile.exists()) {
            if (outputFile.exists()) {
                FileUtils.deleteQuietly(outputFile)
            }
            backupFile.renameTo(outputFile)
        }
    }

    /**
     * Remove filename extension. Removes filename extension from a filename or even full path.
     * If the extension is "gz", it will be removed as well.
     * "/tmp/this" -> "/tmp/this"
     * "/tmp/this.vcf" -> "/tmp/this"
     * "/tmp/this.and.that.vcf" -> "/tmp/this.and.that"
     * "/tmp/this.and.that.vcf.gz" -> "/tmp/this.and.that"
     * "/tmp/this.and.that.vcf.GZ" -> "/tmp/this.and.that"
     * @param filename incoming fulname or even full path
     * @return input without the file extensions (and additionally without .gz if it exists)
     */
    public String filenameWithoutExtension(final String filename) {
        String result = filename
        if (result) {
            while (true) {
                final currentExt = FilenameUtils.getExtension(result).toLowerCase()
                if (!currentExt) {
                    // No extension, we're done
                    break
                }
                result = FilenameUtils.removeExtension(result)
                if (currentExt.toLowerCase() != "gz") {
                    // The removed extension was not "gz". we're done.
                    break
                }
            }
        }
        return result
    }
}

class LuceneExporter extends TsvVcfExporter {

    private static final Logger LOG = Logger.getLogger(LuceneExporter.class);

    final ProgressLogger progress

    Map<String, DbColumnDetails> columnModel
    Document document

    public LuceneExporter() {
        progress = new ProgressLogger(LOG);
        progress.priority = Level.INFO
        progress.displayFreeMemory = true
        progress.itemsName = "rows"
        columnModel = new LinkedHashMap<String, DbColumnDetails>()
        document = new Document()
    }

    /**
     * Import the data into Sqlite.
     * @return 0 - processed correctly
     *          2 - no data to export
     *          5 - exception indexing
     *         99 - unknown status
     */
    @Override
    int convert(final File tsvFile, final VCFParser parser, final String currentTag) {
        //def tsvFile = "UBCIPVC-stats-small.tsv"
        //def filename = "UBCIPVC-stats.tsv.gz"
        //def filename = "c:/temp/UBCIPVC-stats.tsv"

        final String filenameNoExtension = filenameWithoutExtension(tsvFile.getCanonicalPath())
        // remove the tag and any special characters from the table name:
        final String tableName = FilenameUtils.getBaseName(filenameNoExtension).
                replaceFirst(currentTag + "_", "").
                replaceFirst(currentTag + "-", "").replaceAll("[.\\-,!?*%#@^&()]+", "_")
        // we derive the table name automatically from the input file by removing 'tag'_ and the file extension
        System.out.printf("Importing tag=%s filename=%s to tableName=%s%n", currentTag, filenameNoExtension, tableName)
        outputFilename = filenameNoExtension + ".lucene.index"
        outputFile = new File(outputFilename)
        backupFile = backupExisting()

        int numColumns = parser.countAllFields()
        if (numColumns == 0) {
            println "Number of columns in source file is 0"
            return 2
        }

        FileUtils.deleteDirectory(outputFile)

        observeVcfFields(parser)

        def directory = new SimpleFSDirectory(outputFile)

        IndexWriterConfig indexWriterConfig = new IndexWriterConfig(Version.LUCENE_CURRENT, new SpecialSymbolAnalyzer(Version.LUCENE_CURRENT));
        // We are actually using CREATE as we wipe out the directory (above) but we want to use CREATE_OR_APPEND
        // because we've written the metadata file to the directory
        indexWriterConfig.setOpenMode(IndexWriterConfig.OpenMode.CREATE_OR_APPEND)
        IndexWriter writer = new IndexWriter(directory, indexWriterConfig)

        int lineNo = 0
        progress.start("Starting to import ")
        if (DEBUG) {
            println "numColumns=${numColumns} columnModel=${columnModel}"
        }

        int primaryKeyIndex = 0;
        while (parser.hasNextDataLine()) {
            try {
                for (Map.Entry<String, DbColumnDetails> columnDetailsEntry in columnModel.entrySet()) {
                    final DbColumnDetails columnDetails = columnDetailsEntry.value
                    final Fieldable curField = columnDetails.field
                    if (columnDetails.sourceColumnNum == -1) {
                        // -1 is just for primary key. All other values come from database
                        ((NumericField) curField).setIntValue(primaryKeyIndex)
                        columnDetails.minValue = Math.min(columnDetails.minValue, primaryKeyIndex)
                        columnDetails.maxValue = Math.max(columnDetails.maxValue, primaryKeyIndex)
                        primaryKeyIndex++
                    } else {
                        final String curFieldType = columnDetails.columnType
                        final String cellValue = parser.getStringFieldValue(columnDetails.sourceColumnNum) ?: ""
                        final String cellValueLC = cellValue.toLowerCase()
                        if (curFieldType == "Float") {
                            if (!cellValue.isEmpty()) {
                                if (cellValueLC.startsWith("inf") || cellValueLC.startsWith("+inf")) {
                                    ((NumericField) curField).setFloatValue(Float.POSITIVE_INFINITY)
                                    columnDetails.hasPosInf = true
                                } else if (cellValueLC.startsWith("-inf")) {
                                    ((NumericField) curField).setFloatValue(Float.NEGATIVE_INFINITY)
                                    columnDetails.hasNegInf = true
                                } else if (cellValueLC == "nan" || cellValueLC == ".") {
                                    ((NumericField) curField).setFloatValue(Float.NaN)
                                    columnDetails.hasNan = true
                                } else {
                                    final float numericCell = cellValue.toFloat()
                                    ((NumericField) curField).setFloatValue(numericCell)
                                    columnDetails.observeMinMax(numericCell)
                                }
                            } else {
                                // Lucene cannot store "null" or "empty" for a numeric value, so we store
                                // the "special" value of -MIN_VALUE. The app will have to handle this.
                                ((NumericField) curField).setFloatValue(-Float.MIN_VALUE)
                            }
                        } else if (curFieldType == "Integer") {
                            if (!cellValue.isEmpty()) {
                                final int numericCell = cellValue.toInteger()
                                ((NumericField) curField).setIntValue(numericCell)
                                columnDetails.observeMinMax(numericCell)
                            } else {
                                // Lucene cannot store "null" or "empty" for a numeric value, so we store
                                // the "special" value of -MIN_VALUE. The app will have to handle this.
                                ((NumericField) curField).setIntValue(-Integer.MIN_VALUE)
                            }
                        } else {
                            // String or StringSortable
                            curField.setValue(cellValue)
                        }
                    }
                }
                writer.addDocument(document)
            } catch (Exception e) {
                print "  fields="
                dumpList(columnModel.values()*.field)
                e.printStackTrace()
                return 5
            } finally {
                lineNo++
                parser.next()
                if (numberRowsToExport > 0 && lineNo >= numberRowsToExport) {
                    break
                }
                progress.lightUpdate()
            }
        }
        writer.close()
        progress.done()
        println "Indexed ${lineNo} lines"

        // Dump the metadata file first thing
        dumpMetaData()

        return 0
    }

    def dumpMetaData() {
        if (!outputFile.exists()) {
            outputFile.mkdirs()
        }
        def metadataFile = new File(outputFile, "index.metadata.txt")
        println "Writing metadata file ${metadataFile}"
        BufferedWriter writer = metadataFile.newWriter()
        writer.writeLine(["i", "fieldName", "actualFieldName", "fieldType",
                "fieldGroup", "fieldMinValue", "fieldMaxValue",
                "hasPosInf", "hasNegInf", "hasNan"].join("\t"))
        int i = 0
        for (Map.Entry<String, DbColumnDetails> columnDetailsEntry in columnModel.entrySet()) {
            final DbColumnDetails columnDetails = columnDetailsEntry.value

            final String fieldName = columnDetails.dbColumnName
            final String actualFieldName = columnDetails.columnName
            final String fieldType = columnDetails.columnType
            final String fieldGroup = columnDetails.columnGroups.join(",")
            def fieldMinValue = columnDetails.minValue
            def fieldMaxValue = columnDetails.maxValue
            final def hasPosInfValue = columnDetails.hasPosInf
            final def hasNegInfValue = columnDetails.hasNegInf
            final def hasNanValue = columnDetails.hasNan

            if ((fieldType == "Float" && fieldMinValue == Float.MAX_VALUE && fieldMaxValue == Float.MIN_VALUE) ||
                    (fieldType == "Integer" && fieldMinValue == Integer.MAX_VALUE && fieldMaxValue == Integer.MIN_VALUE)) {
                fieldMinValue = 0
                fieldMaxValue = 0
            }

            writer.writeLine([i, fieldName, actualFieldName, fieldType,
                    fieldGroup, fieldMinValue, fieldMaxValue,
                    hasPosInfValue, hasNegInfValue, hasNanValue].join("\t"))
            i++
        }
        writer.close()
        println "...done"
    }

    def observeVcfFields(VCFParser parser) {
        columnModel.clear()
        int numFields = parser.countAllFields()

        final DbColumnDetails pkColumnDetails = new DbColumnDetails()
        columnModel["my_primary_key"] = pkColumnDetails
        pkColumnDetails.columnGroups << "primaryKey"
        pkColumnDetails.columnType = "Integer"
        pkColumnDetails.dbColumnName = "my_primary_key"
        pkColumnDetails.columnName = "my_primary_key"
        pkColumnDetails.columnNum = columnModel.size() - 1
        pkColumnDetails.indexed = true
        pkColumnDetails.sortable = true
        pkColumnDetails.sourceColumnNum = -1  // states this is the primary key
        pkColumnDetails.field = new NumericField(pkColumnDetails.dbColumnName, Field.Store.YES, true)
        pkColumnDetails.defaultMinMaxValues(ColumnType.Integer)

        GroupAssociations groupAssociations = null
        try {
            groupAssociations = parser.getGroupAssociations();
        } catch (NullPointerException e) {
            // It appears TSV may not support group associations?
        }
        for (int i = 0; i < numFields; i++) {
            final String dbColName = "col_${i}"
            final ColumnType fieldType = parser.getFieldType(i)
            final String fieldName = parser.getFieldName(i)

            final DbColumnDetails columnDetails = new DbColumnDetails()
            columnModel[dbColName] = columnDetails
            if (groupAssociations != null) {
                columnDetails.columnGroups.addAll groupAssociations.listGroups(fieldName)
            }
            columnDetails.columnType = fieldType.toString()
            columnDetails.dbColumnName = dbColName
            columnDetails.columnName = fieldName
            columnDetails.columnNum = columnModel.size() - 1
            columnDetails.indexed = true
            columnDetails.sortable = (fieldType != ColumnType.String)
            columnDetails.sourceColumnNum = i
            columnDetails.defaultMinMaxValues(fieldType)
            if (fieldType == ColumnType.Integer || fieldType == ColumnType.Float) {
                columnDetails.field = new NumericField(dbColName, Field.Store.YES, true)
            } else {
                columnDetails.field = new Field(dbColName, "", Field.Store.YES, Field.Index.ANALYZED_NO_NORMS);
            }

            if (fieldType == ColumnType.String) {
                // The sortable version for String columns
                final String sortDbColName = "${dbColName}_s"
                final DbColumnDetails sortColumnDetails = new DbColumnDetails()
                columnModel[sortDbColName] = sortColumnDetails
                sortColumnDetails.columnType = "${fieldType.toString()}Sortable"
                sortColumnDetails.dbColumnName = sortDbColName
                sortColumnDetails.columnName = sortDbColName
                sortColumnDetails.columnNum = columnModel.size() - 1
                sortColumnDetails.indexed = true
                sortColumnDetails.sortable = true
                sortColumnDetails.field = new Field(sortDbColName, "", Field.Store.NO, Field.Index.NOT_ANALYZED_NO_NORMS);
                sortColumnDetails.defaultMinMaxValues(ColumnType.String)
                sortColumnDetails.sourceColumnNum = i
            }
        }

        // Add all of the fields we've created to the document.
        for (final Map.Entry<String, DbColumnDetails> entry in columnModel) {
            document.add(entry.value.field)
            if (DEBUG) {
                println "Adding field to document ${entry.value.field}"
            }
        }
    }
}

class SqliteExporter extends TsvVcfExporter {

    private static final Logger LOG = Logger.getLogger(SqliteExporter.class);

    final static DB_DOUBLE_TYPE = "BINARY_FLOAT"  //BINARY_DOUBLE for oracle
    final static DB_STRING_TYPE = "VARCHAR2(2)" //VARCHAR2(50) for oracle, just VARCHAR for hsqldb
    final static String DB_INTEGER_TYPE = "PLS_INTEGER"  // Signed Integer for oracle

    final ProgressLogger progress
    final Sqlite sqliteDatabaseService
    final StringBuilder sql

    Sql db = null

    public SqliteExporter() {
        sqliteDatabaseService = new Sqlite()
        sql = new StringBuilder()
        progress = new ProgressLogger(LOG);
        progress.displayFreeMemory = true;
        progress.priority = Level.INFO
        progress.itemsName = "rows"
    }

    void openSqliteDatabase() {
        println "Obtaining sql connection"
        db = sqliteDatabaseService.borrowConnection(outputFilename)
        // We don't need transactions since we write once and query many times, disable the journal.
        db.connection.prepareStatement("PRAGMA main.page_size = 4096;").execute();
        db.connection.prepareStatement("PRAGMA cache_size = 100000;").execute();
        db.connection.prepareStatement("PRAGMA main.locking_mode=EXCLUSIVE;").execute();
        db.connection.prepareStatement("PRAGMA synchronous = OFF;").execute();
        db.connection.prepareStatement("PRAGMA journal_mode = OFF;").execute();
        db.connection.prepareStatement("PRAGMA main.temp_store = MEMORY;").execute();
    }

    /**
     * Import the data into Sqlite.
     * @return 0 - processed correctly
     *          1 - import file not found
     *          2 - No columns (fields) found in input file
     *          3 - SQLException importing data
     *          4 - Other Exception importing data
     *         99 - unknown status
     */
    public int convert(final File tsvFile, final VCFParser parser, final String currentTag) {
        int retval = 99 // Unknown result
        Boolean oldAutoCommit = null

        println "Importing data from ${tsvFile} to Sqlite"
        exec.queueMessage currentTag, "Importing data into database"
        int position
        int tabPosition
        int numLinesRead = 0
        int goodDataLines = 0
        boolean lineComplete
        String cellValue
        String cellValueLC
        int startPosition
        try {
            int numColumns = parser.countAllFields()
            if (numColumns == 0) {
                println "Number of columns in source file is 0"
                retval = 2
            }

            final String filenameNoExtension = filenameWithoutExtension(tsvFile.getCanonicalPath())
            // remove the tag and any special characters from the table name:
            final String tableName = FilenameUtils.getBaseName(filenameNoExtension).
                    replaceFirst(currentTag + "_", "").
                    replaceFirst(currentTag + "-", "").replaceAll("[.\\-,!?*%#@^&()]+", "_")
            // we derive the table name automatically from the input file by removing 'tag'_ and the file extension
            System.out.printf("Importing tag=%s filename=%s to tableName=%s%n", currentTag, filenameNoExtension, tableName)
            outputFilename = filenameNoExtension + ".db"
            outputFile = new File(outputFilename)
            backupFile = backupExisting()

            openSqliteDatabase()

            // If the table exists already, drop it, we'll recreate.
            sqliteDatabaseService.dropTable db, tableName
            sqliteDatabaseService.dropTable db, "${tableName}_colnames"

            oldAutoCommit = db.connection.autoCommit
            println "Setting autoCommit to false"
            db.connection.autoCommit = false
            println "Importing data into ${tableName}"

            def columns = createTableSchema(db, parser, tableName);

            // Create the indexes first, to try and avoid multiple passes through the data (one for each index):
            boolean createIndexes = true

            progress.start("Starting to import ")
            final StringBuilder insertSql = new StringBuilder()

            insertSql << "INSERT INTO " << tableName << " values ("
            for (int globalFieldIndex = 0; globalFieldIndex < numColumns + 1; globalFieldIndex++) {
                if (globalFieldIndex != 0) {
                    insertSql << ","
                }
                // value to be taken from row:
                insertSql << "?"
            }
            insertSql << ")"
            System.out.println("insertSql: " + insertSql);

            PreparedStatement ps = db.connection.prepareStatement(insertSql.toString())
            int primaryKeyIndex = 0;
            while (parser.hasNextDataLine()) {

                ps.setInt(1, primaryKeyIndex++);
                for (int globalFieldIndex = 0; globalFieldIndex < numColumns; globalFieldIndex++) {

                    int index = globalFieldIndex + 2;
                    cellValue = parser.getStringFieldValue(globalFieldIndex)
                    cellValueLC = cellValue.toLowerCase()
                    if ("".equals(cellValue)) {
                        // missing value.
                        ps.setString(index, "")
                    } else {
                        if (parser.getFieldNumValues(globalFieldIndex) == 1) {

                            switch (parser.getFieldType(globalFieldIndex)) {

                                case ColumnType.Float:
                                    if (cellValueLC.indexOf("inf") == 0) {
                                        ps.setFloat(index, Float.POSITIVE_INFINITY)
                                    } else if (cellValueLC.indexOf("+inf") == 0) {
                                        ps.setFloat(index, Float.POSITIVE_INFINITY)
                                    } else if (cellValueLC.indexOf("-inf") == 0) {
                                        ps.setFloat(index, Float.NEGATIVE_INFINITY)
                                    } else if (cellValueLC.indexOf("nan") == 0) {
                                        ps.setFloat(index, Float.NaN)
                                    } else if (cellValueLC.equals(".")) {
                                        ps.setFloat(index, Float.NaN)
                                    } else {
                                        ps.setFloat(index, cellValue.toFloat())
                                    }
                                    break;

                                case ColumnType.Integer:

                                    ps.setInt(index, cellValue.toInteger())
                                    break;
                                default:
                                    ps.setString(index, cellValue)
                                    break;
                            }
                        } else {
                            // multiple values are stored as String
                            ps.setString(index, cellValue)
                        }
                    }
                }
                ps.addBatch()

                goodDataLines++
                numLinesRead++

                if (numLinesRead % 10000 == 0) {
                    ps.executeBatch()
                    db.commit()
                    // println "... ${String.format("%,d", numLinesRead)} rows comitted"

                    if (numberRowsToExport > 0 && goodDataLines >= numberRowsToExport) {
                        break
                    }
                }

                parser.next()
                progress.lightUpdate()
            }
            ps.executeBatch()
            db.commit()
            progress.done()
            println "... Inserts completed, a total of " << String.format("%,d", numLinesRead) << " rows"
            db.commit()
            println "... Done"

            if (createIndexes) {
                progress.expectedUpdates = columns.size()
                progress.start("Starting to index")
                println "Indexing.."

                def sqlIndex = new StringBuilder()
                sqlIndex << "CREATE INDEX " << tableName << "_PRIMARY ON " << tableName << " (my_primary_key);"

                for (int i = 0; i < columns.size(); i++) {
                    sqlIndex.setLength(0)
                    ColumnType type = parser.getFieldType(i);
                    boolean indexColumn = false;
                    switch (type) {
                        case ColumnType.Float:
                        case ColumnType.Character:
                            indexColumn = true;
                    }
                    if (i < 5) {
                        // always index the first 5 columns
                        indexColumn = true
                    }
                    if (indexColumn) {
                        sqlIndex << "CREATE INDEX " << tableName << "_" << i << " ON " << tableName << " (" << columns[i] << ");"
                    }
                    System.out.println("indexing: " + sqlIndex);
                    db.connection.createStatement().executeUpdate(sqlIndex.toString())
                    db.connection.commit()
                    progress.lightUpdate()
                }
                println "Done.."
                progress.done()
            }

            retval = 0
        } catch (SQLException e) {
            LOG.error("SQLException importing diffexp data into database", e)
            retval = 3
            exec.queueMessage currentTag, "Error importing data into database"
            throw e
        } catch (Exception e) {
            LOG.error("Other exception importing diffexp data into database", e)
            retval = 4
            exec.queueMessage currentTag, "Error importing data into database"
            throw e
        } finally {
            if (db != null) {
                db.commit()
                if (oldAutoCommit != null) {
                    db.connection.autoCommit = oldAutoCommit
                }
                sqliteDatabaseService.returnConnection(db)
                if (retval == 0) {
                    exec.queueMessage currentTag, "Data imported into database"
                } else {
                    restoreBackup()
                }
            }
        }

        return retval
    }

    private Object createTableSchema(Sql db, VCFParser parser, String tableName) {
        def sql = new StringBuilder()
        def insertSql = new StringBuilder()
        println "Creating tables"
        def colTypes = []
        def columnNames = []
        def colGroups = []
        def columns = []
        int numFields = parser.countAllFields()
        GroupAssociations groupAssociations = parser.getGroupAssociations();
        sql << "CREATE TABLE " << tableName << " (my_primary_key Integer PRIMARY KEY, "
        insertSql << "INSERT INTO " << tableName << " values ("
        for (int i = 0; i < numFields; i++) {
            if (i != 0) {
                insertSql << ","
                sql << ","
            }
            String colName = parser.getFieldName(i);

            columns << "col_${i}"
            columnNames << colName
            colGroups << groupAssociations.listGroupsAsString(colName)
            if (parser.getFieldNumValues(i) == 1) {

                switch (parser.getFieldType(i)) {

                    case ColumnType.Float:
                        sql << "${columns[i]} ${DB_DOUBLE_TYPE}"   // NUMBER, BINARY_DOUBLE
                        colTypes << "float"
                        break;

                    case ColumnType.Integer:
                        sql << "${columns[i]} ${DB_INTEGER_TYPE}"   // NUMBER, BINARY_DOUBLE
                        colTypes << "int"
                        break;
                    default:
                        sql << "${columns[i]} ${DB_STRING_TYPE}"
                        colTypes << "string"
                        break;

                }
            } else {
                // multiple values are stored as string
                sql << "${columns[i]} ${DB_STRING_TYPE}"
                colTypes << "string"
            }
            insertSql << "?"
        }
        sql << ")"
        insertSql << ")"
        insertSql = insertSql.toString()
        println sql.toString()
        println insertSql
        db.execute sql.toString()
        println "... storage table created"

        sql.setLength(0)
        sql << "CREATE TABLE ${tableName}_colnames (dbColName ${DB_STRING_TYPE}, screenColName ${DB_STRING_TYPE}, colType ${DB_STRING_TYPE}, colGroups ${DB_STRING_TYPE})"
        println sql.toString()
        db.execute sql.toString()
        sql.setLength(0)
        sql << "INSERT INTO ${tableName}_colnames values ('my_primary_key', 'row', 'integer', 'primary_key')"
        println sql.toString()
        db.execute sql.toString()
        for (int i = 0; i < numFields; i++) {
            sql.setLength(0)
            sql << "INSERT INTO " << tableName << "_colnames values ('${columns[i]}', '${columnNames[i]}', '${colTypes[i]}', '${colGroups[i]}')"
            println sql.toString()
            db.execute sql.toString()
        }
        db.commit()
        println "... table of column names created"
        return columns
    }
}
