#!/usr/bin/env groovy

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

//input args are filenames to combine

class FileData{
    Reader input

    String line
    String[] parts
}

def fileList = args.collect {
    def out = new FileData()
    out.input = new File(it).newReader()
    out.input.readLine()
    out
}

def read = {
    fileList.each {
        it.line = it.line ?: it.input.readLine()
        it.parts = it.line?.split("\t")
    }
    fileList
}

def window

while ((window = read()).any {it.line != null}){
    def minPosition = window.findAll {it.line != null}.min {
        it.parts[1].toInteger()
    }.parts[1].toInteger()

    def chromosome = ""

    def row = window.collect {
        if(it.line != null && it.parts[1].toInteger() == minPosition){
            it.line = null
            chromosome = it.parts[0]
            it.parts[2..4].join("\t")
        }
        else{
            ["NA", "NA", "NA"].join("\t")
        }
    }

    println "${chromosome}\t${minPosition}\t${row.join("\t")}"

}