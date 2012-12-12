/**
 * Script that extracts counts from the STAR splicing output, and formats the counts in such a way that they can be fed
 * to EdgeR to estimate p-values of differential splicing usage.
 * @author Fabien Campagne
 *         Date: Jul / 25 / 12
 *         Time: 12: 50 PM
 */

import collection.mutable.ListBuffer
import it.unimi.dsi.fastutil.ints.IntArrayList
import it.unimi.dsi.io.{FastBufferedReader, LineIterator}
import it.unimi.dsi.io.{LineIterator, FastBufferedReader}
import java.io.FileReader
var nextSampleIndex = -1;
val sampleMap = scala.collection.mutable.Map[String, Int]()
val indexToSample = scala.collection.mutable.Map[Int, String]()
for (file <- args) {
  val lines: LineIterator = new LineIterator(new FastBufferedReader(new FileReader(file)))
  //val lines: Iterator[String] = scala.io.Source.fromFile(file).getLines()
  // ignore the header line:
  header = lines.next().toString
  while (lines.hasNext) {
    val line = lines.next()
    val tokens = line.toString.split("\t")

    val sample = tokens(0)
    val sampleIndex :Int = if (sampleMap.contains(sample)) {
      sampleMap.get(sample).get
    } else {
      nextSampleIndex+=1
      val newSampleIndex = (nextSampleIndex)
      sampleMap += (sample -> newSampleIndex)
      newSampleIndex
    }
    sampleMap += (sample -> sampleIndex)
    indexToSample += (sampleIndex -> sample)
  }
}

//printf("Found %d samples%n",sampleMap.size)
//println("Executing with arguments: " + args.mkString(","))
// data map: key is site, at each site, an array with a count for each sample (by sampleIndex)
val dataMap = scala.collection.mutable.Map[String, IntArrayList]()

var header: String = ""
val COUNT_COL_INDEX = 7

for (file <- args) {
  val lines: LineIterator = new LineIterator(new FastBufferedReader(new FileReader(file)))
  //val lines: Iterator[String] = scala.io.Source.fromFile(file).getLines()
  // ignore the header line:
  header = lines.next().toString
  while (lines.hasNext) {
    val line = lines.next()
    var site: String = ""
    val tokens = line.toString.split("\t")

    val sample:String = tokens(0)
    val sampleIndex: Int=sampleMap.getOrElse(sample,-1)
    site= tokens(1)+":" +tokens(2)+"-"+tokens(3)+":"+tokens(4)

    val spliceJunctionCount = Integer.parseInt(tokens(6))
    val newList:IntArrayList=new IntArrayList()
    newList.size(sampleMap.size)
    val dataForSamplesAtSite: IntArrayList=dataMap.getOrElse(site, newList)
//   printf("Processing sampleIndex=%d sample=%s%n",sampleIndex,sample)
    dataForSamplesAtSite.set(sampleIndex,spliceJunctionCount.toInt)
    dataMap += (site -> dataForSamplesAtSite)
//    printf("Processing sampleIndex=%d sample=%s list=%s %n",sampleIndex,sample, dataForSamplesAtSite)
  }
}

//printf("chromosome\tintron-first-base\tintron-last-base\t%s",sampleMap.keys.mkString("\t"))
//printf("chromosome|intron-first-base|intron-last-base|\telement-type\t%s",sampleMap.keys.mkString("\t"))
val sortedKeys= ListBuffer[String]()
for (sampleIndex <- 0 to sampleMap.size-1) {
  sortedKeys += indexToSample.get(sampleIndex).get
}
printf("element-id\telement-type\t%s",sortedKeys.mkString("\t"))
println()
val LOG2: Double = StrictMath.log(2)
var i = 0;
for (site: String <- dataMap.keys) {

  val siteTokens: Array[String] = site.split("[|]")

 // print(siteTokens.mkString("\t"))
  print(site)
  print ("\tSPLICE\t")
  val value=dataMap.get(site)
 // println("value: "+value)
  val siteMap: IntArrayList=value match {
    case Some(l) => l
    case None =>new IntArrayList()
  }
//  println("siteMap: "+siteMap)
  for (sampleIndex <- 0 to sampleMap.size-1) {
      print(siteMap.getInt(sampleIndex))
      if (sampleIndex!=sampleMap.size-1) print("\t")
  }
  println()
  /*
  } */

  //if (i>2) System.exit(1)
  //i+=1
}

