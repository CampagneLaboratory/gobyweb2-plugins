/*
 * Copyright (c) 2011-2012  by Cornell University  and the  Cornell Research Foundation,
 * Inc. All Rights Reserved.
 *
 * GobyWeb plugins  are released  under the LGPL3 license,  unless  otherwise specified
 * by the license of a specific plugin. See the file LGPL3.license in this distribution
 * for a copy of the LGPL license.
 *
 * When a plugin is not released under the LGPL3 license,  the comments  at the top  of
 * the plugin's config.xml will indicate how that specific plugin is released/licensed.
 */

/**
 * Make an output file which maps alignments to reads.
 * @param gobywebObj the gobyweb object just configured for this plugin
 * @param tempDir temporary directory where to write files that will be transferred to the cluster with the plugin
 * @return exit code, 0 means executed normally
 */
int execute(final Object gobywebObj, final File tempDir, final Map bindings) {
    final File outputFile = new File(tempDir, "splicejunctioncoverage-groups.tsv")
    final File sampleGroupsFile = new File(tempDir, "sampleGroups.tsv")
    final PrintWriter writer = outputFile.newPrintWriter()
    final PrintWriter sampleGroupsWriter = sampleGroupsFile.newPrintWriter()
    Object pathService = bindings.get("pathService")
    Object config = bindings.get("config")

    try {
        Map<String, String> grpToName = gobywebObj.grpToName

        gobywebObj.grpToAligns.each { String key, Object alignment ->
            println "alignment=" + alignment
            println "alignJob=" + alignment.alignJob
            String credentials = config.gobyweb.webServerSshPrefix
            String path = pathService.usersExistingWebJobResultsDir(alignment.alignJob)
            final String compactReadsLocation = (alignment.alignJob.sample.compactReads as List)[0].url.split(":")[1]

            // the key has the form groupIndex - numAlignmentsInGroup. We extract the group index:
            String groupIndex = key.split("-")[0]
            // we find the group name in the grpToName map:
            final String groupName = grpToName.get(groupIndex)
            // we get the alignment basename:
            final String spliceFilename = "${alignment.alignJob.tag}-SpliceJunctionCoverage-all.tsv"
            // we write the tab delimited output:
            writer.println "${credentials}:${path}/${spliceFilename}\t${groupName}"
            sampleGroupsWriter.println "${alignmentBasename(alignment)}\t${groupName}"
        }
    } finally {
        writer.close()
        sampleGroupsWriter.close()
    }
    return 0
}

/**
 * This comes from alignmentService.alignmentFilename().
 * @param alignment the alignment in question
 * @param extension the extension to provide a filename for
 * @return
 */
public String alignmentFilename(Object alignment, String extension) {
    // some versions of GobyWeb stored "tag-basename" in the basename.
    if (alignment.basename.startsWith(alignment.alignJob.tag)) {
        return "${alignment.basename}.${extension}"
    } else {
        // if the basename does not include the tag, make sure we return a filename that includes it:
        return "${alignment.alignJob.tag}-${alignment.basename}.${extension}"
    }
}


/**
 * This comes from alignmentService.alignmentFilename().
 * @param alignment the alignment in question
 * @return basename of the alignment.
 */
public String alignmentBasename(Object alignment) {

    // some versions of GobyWeb stored "tag-basename" in the basename.
    if (alignment.basename.startsWith(alignment.alignJob.tag)) {
        return "${alignment.basename}"
    } else {
        // if the basename does not include the tag, make sure we return a filename that includes it:
        return "${alignment.alignJob.tag}-${alignment.basename}}"
    }
}
