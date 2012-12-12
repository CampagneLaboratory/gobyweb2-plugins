/**
 * Script that adds gene ids to a junction output file. The script extracts the intron boundaries from the first
 * column of the input, determines which gene fully contains these boundaries and adds the gene id as first column
 * of the output.
 *
 *
 * @author Fabien Campagne
 *         Date: Jul / 27 / 12
 *         Time: 12: 50 PM
 */

import it.unimi.dsi.fastutil.ints.IntArrayList
import it.unimi.dsi.io.{FastBufferedReader, LineIterator}
import it.unimi.dsi.io.{LineIterator, FastBufferedReader}
import java.io.FileReader
import edu.cornell.med.icb.goby.algorithmic.data.Interval

val annotations = new edu.cornell.med.icb.goby.algorithmic.algorithm.RandomAccessAnnotations()
val annotationFilename = args(0)
//System.out.println("Loading annotations:" + annotationFilename)
//System.out.flush()
annotations.loadAnnotations(annotationFilename)
//System.out.println("Done loading annotations:" + annotationFilename)
val emptyInterval = new Interval
emptyInterval.id = ""

for (file <- args.tail) {
  val lines: LineIterator = new LineIterator(new FastBufferedReader(new FileReader(file)))

  val header = lines.next().toString
  val firstToken=header.split("\t")(0)
  System.out.println(header.replaceFirst(firstToken, "gene-id\tchromosome\tintronFirstBase\tintronLastBase\tstrand"))
  while (lines.hasNext) {
    val line = lines.next()
    val tokens = line.toString.split("\t")
    val ids = tokens(0).split("[:]")
    val chromosome = ids(0)
    val start_end=ids(1).split("[-]")
    val strand=ids(2)
    val start = Integer.parseInt(start_end(0))
    val end = Integer.parseInt(start_end(1))
    var interval:Interval = annotations.find(chromosome, start, end)
    if (interval == null) {
      interval = emptyInterval
    }
    val newIds=List(chromosome,start,end,strand)
    System.out.printf("%s\t%s\t%s%n", interval.id, newIds.mkString("\t"), line.substring(line.indexOf('\t')+1))
  }
  System.out.flush()

}

