#! /usr/bin/env python

# Copyright 2010, 2011 Martin C. Frith

import fileinput, itertools, optparse, os, signal, sys

def batches(lines):
    for line in lines:
        if line.startswith("# batch"):
            yield line
        else:
            print line,
    while True:
        yield None

def lastMergeBatches(fileNames):
    files = map(fileinput.input, fileNames)
    b = map(batches, files)

    for i in itertools.izip(*b):
        j = filter(None, i)
        if j: print j[0],
        else: break

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = "%prog files"

    description = "Read files of lastal output, merge corresponding batches, and write them."

    op = optparse.OptionParser(usage=usage, description=description)
    opts, args = op.parse_args()
    if not args: op.error("please give me some file names")

    try: lastMergeBatches(args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
