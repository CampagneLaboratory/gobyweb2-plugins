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

//aggregate methylation calls across positions
//use Red-Black Tree Map built in java implementation
//arg 1 = path to methylation calls

println args

def map
def activeChromosome

def MAX_RANGE = 10000
def TRIM_RANGE = 1000

def exportAll = {
    map.each { k, v ->
        println "${activeChromosome}\t${k}\t${v[0]}\t${v[1]}\t${v[0].doubleValue() / (v[0] + v[1])}"
    }
}

new File(args[0]).splitEachLine("\t") { fields ->
    if(fields.size() != 5)
        return

    def chromosome = fields[2]

    if(chromosome != activeChromosome){
        exportAll()
        map = new TreeMap<Integer, Integer[]>()
    }
    activeChromosome = chromosome

    def position = fields[3].toInteger()

    def base = map.get(position, [0, 0])

    (fields[1] == '+' ? base[0]++ : base[1]++)

    if(position - map.firstKey() >= MAX_RANGE){

        def export = map.headMap(position - TRIM_RANGE)
        map = map.tailMap(position - TRIM_RANGE)

        export.each {k, v ->
            println "${activeChromosome}\t${k}\t${v[0]}\t${v[1]}\t${v[0].doubleValue() / (v[0] + v[1])}"
        }
    }
}

map.each {k, v ->
    println "${activeChromosome}\t${k}\t${v[0]}\t${v[1]}\t${v[0].doubleValue() / (v[0] + v[1])}"
}