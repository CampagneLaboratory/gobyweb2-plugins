
import org.campagnelab.mercury.cli.JobInterface

//broker host
def hostname = args[0]

//broker port
def port = args[1]

def tag = args[2]

def jobDir = args[3]

def traceMap = [:]

def warningLines = []

// Number of lines to keep in buffer
def BUFFER_SIZE = 100

// Pattern for stack trace line
def TRACE_LINE_PATTERN = '^[\\s\\t]+at .*$'

// Log line pattern between which we try to capture full trace
def LOG_LINE_PATTERN = '^([<#][^/]|\\d\\d).*$'

// List of patterns to replace in final captured stack trace line 
// (e.g. replace date and transaction information that may make similar traces to look as different)
def REPLACE_PATTERNS = [
        '^\\d+-\\d+\\@.*?tksId: [^\\]]+\\]',
        '^<\\w+ \\d+, \\d+ [^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <',
        '^####<[^>]+?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <[^>]*?> <',
        '<([\\w:]+)?TransaktionsID>[^<]+?</([\\w:]+)?TransaktionsID>',
        '<([\\w:]+)?TransaktionsTid>[^<]+?</([\\w:]+)?TransaktionsTid>'
]

def NO_SUCH_FILE_LINE_PATTERN = '.*: No such file or directory$'

new File('.').eachFile { File file ->
    if (file.name.contains('.log') || file.name.contains('.out') || file.name.matches("[A-Z]{7}\\..*")) {

        //  println "Scanning filename="+file.name
        def bufferLines = []

        file.withReader { Reader reader ->
            int lineUpCount = 0
            while (reader.ready()) {
                def String line = reader.readLine()
                if (line.matches(TRACE_LINE_PATTERN)) {
                    //println "Matched trace line pattern: "+line
                    def trace = []
                    for (def i = bufferLines.size() - 1; i >= 0; i--) {
                        lineUpCount++
                        if (!bufferLines[i].matches(LOG_LINE_PATTERN)) {
                            //   println "Did not match LOG_LINE_PATTERN: "+bufferLines[i]

                            trace.add(0, bufferLines[i])
                            if (lineUpCount > 5) {
                                break
                            };
                        } else {
                            //     println "Matched LOG_LINE_PATTERN: "+bufferLines[i]
                            trace.add(0, bufferLines[i])
                            break
                        }
                    }
                    trace.add(line)
                    if (reader.ready()) {
                        line = reader.readLine()
                        while (line.matches(TRACE_LINE_PATTERN)) {
                            trace.add(line)
                            if (reader.ready()) {
                                line = reader.readLine()
                            } else {
                                break;
                            }
                        }
                    }
                    lineUpCount = 0
                    def traceString = trace.join("\n")
                    REPLACE_PATTERNS.each { pattern ->
                        traceString = traceString.replaceAll(pattern, '')
                    }
                    if (traceMap.containsKey(traceString)) {
                        traceMap.put(traceString, traceMap.get(traceString) + 1)
                    } else {
                        traceMap.put(traceString, 1)
                    }
                } else if (line.matches(NO_SUCH_FILE_LINE_PATTERN)) {
                    warningLines.add(line)
                }
                // Keep the buffer of last lines.
                bufferLines.add(line)
                if (bufferLines.size() > BUFFER_SIZE) {
                    bufferLines.remove(0)
                }

            }
        }
    }
}

traceMap = traceMap.sort { it.value }

traceMap.reverseEach { trace, number ->
    def args = [
            "--broker-hostname", hostname,
            "--broker-port", port,
            "--job-tag", tag,
            "--description", trace,
            "--phase", "post_process",
            "--category", "ERROR",
            "--jndi-config", "${jobDir}/mercury.properties"]
    JobInterface.processAPI(args as String[])
}

warningLines.each { line ->
    def args = [
            "--broker-hostname", hostname,
            "--broker-port", port,
            "--job-tag", tag,
            "--description", line,
            "--phase", "post_process",
            "--category", "ERROR",
            "--jndi-config", "${jobDir}/mercury.properties"]
    JobInterface.processAPI(args as String[])
}

System.exit(0)