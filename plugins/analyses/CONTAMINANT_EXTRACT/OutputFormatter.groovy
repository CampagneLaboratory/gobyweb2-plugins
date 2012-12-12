#!/usr/bin/env groovy

def ARGS_LENGTH = 6
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

//arg 1: path to accession-name map file
//arg 2: path to input file
//arg 3: path to full output file
//arg 4: path to summary output file
//arg 5: E-value threshold
//arg 6: Identity threshold

println args

def eValueThreshold = args[4].toDouble()
def identityThreshold = args[5].toDouble()

def printUsage(){
    println "Incorrect Syntax"
    println "./OutputFormatter.groovy map_file input full_out summ_out e_value_thresh ident_thresh"
}

if(args.length != ARGS_LENGTH){
    printUsage()
    System.exit(1)
}

//read in the names file
def nameMap = [:]
new File(args[0]).splitEachLine("\t"){
    nameMap[it[0]] = it[1]
}


def outFull = new File(args[2])
def sampleMap = [:]

outFull.write("Contaminant Species\tAccession Number\tSample\tContig\tAlignment Size\tScore\tPercent Identity\tE-value\n")

new File(args[1]).splitEachLine("\t"){

    def percentIdentity = (it[4].toDouble() / it[3].toDouble()) * 100
    outFull << nameMap[it[0].trim()] << "\t" <<
            it[0..4].join("\t") << "\t" <<
            percentIdentity << "\t" <<
            it[5] << "\n";

    //get existing organism list for sample, or make empty one if it doesnt exist yet
    def orgList = sampleMap.get(it[1], [])

    //add organism to list if it passes filter
    if(it[5].toDouble() < eValueThreshold && percentIdentity > identityThreshold){
        orgList << nameMap[it[0].trim()]
    }
}


def outsumm = new File(args[3])
outsumm.write("Sample\tOrganism\tNum Significant Matches\n")
//this is where I will make the summary output

sampleMap.collect { sample, List orgList ->
    orgList.unique(false).collect { [sample, it,  orgList.count(it)] }  //collect counts for each sample-organism pair
}.inject([]) { acc, ele -> //flatten one level down
    acc.addAll(ele)
    acc
}.each { //write each record to output file
    outsumm << it.join("\t") << "\n"
}