
#These variables are here only for standalone testing
#JOB_DIR=.
#RESOURCES_GROOVY_EXECUTABLE=groovy
#TMPDIR=.

if [ ! -x "${JOB_DIR}" ]; then
       echo "You must define JOB_DIR before running this job."
       exit 1;
fi


if [ ! -x "${TMPDIR}" ]; then
       echo "You must define TMPDIR before running this job."
       exit 1;
fi

. ${JOB_DIR}/constants.sh;
. ${JOB_DIR}/auto-options.sh;
. ${TMPDIR}/exports.sh

export JAVA_OPTS="-Xms1536m -Xmx3072m"
${RESOURCES_GROOVY_EXECUTABLE} ${JOB_DIR}/ArtifactDownloader.groovy "$@"