#! /bin/sh

# This script reads MAF-format alignments with lastal header
# information, removes "uninteresting" alignments, and writes the
# remainder.  Specifically, it removes alignments that are "dominated"
# in both genomes.  See last-remove-dominated.py for the definition of
# "dominated".  This procedure is likely to remove many paralogs, and
# it is unlikely to remove one-to-one orthologs.

# If option "-d" is specified, then it removes alignments that are
# dominated in either genome.

sortOpt=
while getopts hd opt
do
    case $opt in
	h)  cat <<EOF
Usage: $(basename $0) [options] last-output.maf

Options:
  -h  show this help message and exit
  -d  reduce alignments more aggressively
EOF
	    exit
	    ;;
	d)  sortOpt="-d"
	    ;;
    esac
done
shift $((OPTIND - 1))

PATH=$PATH:$(dirname $0)  # assume the other scripts are in the same directory

{
    # remove alignments that are dominated in the upper sequence
    maf-sort.sh "$@" |
    last-remove-dominated.py

    # remove alignments that are dominated in the lower sequence
    maf-swap.py "$@" |
    maf-sort.sh |
    last-remove-dominated.py |
    maf-swap.py

} | maf-sort.sh $sortOpt  # merge the results
