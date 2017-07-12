#!/bin/bash

. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh

EVENT_FILE=${TMPDIR}/events-`date +%s`.proto
java ${JAVA_OPTIONS} -Dlog4j.configuration=${SGE_O_WORKDIR}/log4j.properties -jar ${SGE_O_WORKDIR}/events.jar \
       org.campagnelab.gobyweb.events.tools.AppendEvent -p ${EVENT_FILE} --message "$@"


if [ -z "${QUEUE_WRITER_POSTFIX}" ]; then

  java ${JAVA_OPTIONS} -Dlog4j.configuration=${SGE_O_WORKDIR}/log4j.properties -jar ${SGE_O_WORKDIR}/events.jar \
       org.campagnelab.gobyweb.events.tools.PushEvents -p ${EVENT_FILE} ${QUEUE_WRITER_POSTFIX} "$@"

else
  echo "Status update disabled: $*"
fi

