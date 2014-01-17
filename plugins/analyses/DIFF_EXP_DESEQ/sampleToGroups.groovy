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
 * Make an output file which maps sample name to groups. The file contains two columns, separated by a tab character.
 * The first column is the sample name, which is also the basename of the alignment file for this sample.
 * The second column is the group name.
 * @param gobywebObj the gobyweb object just configured for this plugin
 * @param tempDir temporary directory where to write files that will be transferred to the cluster with the plugin
 * @return exit code, 0 means executed normally
 */
int execute(final Object gobywebObj, final File tempDir) {
    final File outputFile = new File(tempDir, "sampleToGroups.tsv")
    final PrintWriter writer = outputFile.newPrintWriter()
    try {

        Map<String, String> grpToName=gobywebObj.grpToName
        gobywebObj.grpToAligns.each { String key, Object alignment ->
            // the key has the form groupIndex - numAlignmentsInGroup. We extract the group index:
            String groupIndex=key.split("-")[0]
            // we find the group name in the grpToName map:
            final String groupName=grpToName.get(groupIndex)
            // we get the alignment basename:
            final String alignJobBasename = "${alignmentBasename(alignment)}"
            // we write the tab delimited output:
            writer.println "${alignJobBasename}\t${groupName}"
        }
    } finally {
        writer.close()
    }
    return 0
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
