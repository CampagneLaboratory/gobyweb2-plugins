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

/**
 * Updates NGA of the status of a job based on the command line.
 */

class QueueWriter {

    String queueMessageDir

    public static void main(final String[] args) {
       new QueueWriter().run(args)
    }

    def run(final String[] args) {
        def opt = parseCommandLine(args)
        if (!opt) {
            return
        }
        def queueItem = new QueueWriterQueueItem(opt)

        File messageFile = File.createTempFile("${queueItem.creationDate}-msg-${queueItem.handlerService}",".queue-message")
        messageFile.write queueItem.toString()
        String command = "scp ${messageFile} ${queueMessageDir}/"
        command.execute().waitFor()
        messageFile.delete()
    }

    def parseCommandLine(final String[] args) {
        def cli = new CliBuilder(usage:'QueueWriter.groovy')
        cli.t(longOpt:'tag', required:true, args:1, "tag of QueueItem")
        cli.s(longOpt:'status', required:true, args:1, "start, failed, completed depending on status of task or job")
        cli.d(longOpt:'description', required:true, args:1, "description of QueueItem")
        cli.i(longOpt:'index', required:true, args:1, "The task index, 1 based, includes all sub tasks")
        cli.j(longOpt:'job-type', required:true, args:1, "The job type (index, align, concat, job)")
        cli.q(longOpt:'queue-message-dir', required:true, args:1,
            "An alternative tp http-prefix, specifies where to store queue message. " +
            "Generally this is somethihng like username@servername:directory")
        cli.r(longOpt:'handler-service', required:false, args:1, 
            "The service that will handle the message (sampleService, diffExpService, alignJobService)")
        def opt = [:]
        def options = cli.parse(args)
        if (options) {
            opt.status = options.s
            opt.description = options.d
            opt.tag = options.t
            opt.index = stringToInteger(options.i)
            opt.jobType = options.j
            queueMessageDir = options.q
            opt.handlerService = options.r
        }
        return opt
    }

    def stringToBoolean(str) {
        if ((str == null) || (str.length() == 0)) {
            return false
        }
        def firstChar = (str.toLowerCase())[0]
        return firstChar == "y" || firstChar == "t"
    }

    def stringToInteger(str) {
        if ((str == null) || (str.length() == 0)) {
            return 0
        }
        if (str.isInteger()) {
            return str as Integer
        } else {
            return 0
        }
    }
}

public class QueueWriterQueueItem {
    long creationDate = (new Date()).time
    String tag
    String status
    String description
    int index
    String jobType
    String handlerService
    public String toString() {
        "[ tag:'${tag}', creationDate:${creationDate}, status:'${status}', " +
        "description:'${description}', index:${index}, jobType:'${jobType}', " +
        "handlerService:'${handlerService}' ]"
    }
}