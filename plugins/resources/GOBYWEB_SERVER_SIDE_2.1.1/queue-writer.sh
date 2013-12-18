#!/bin/bash

. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh

if [ -z "${QUEUE_WRITER_POSTFIX}" ]; then
  ${SGE_O_WORKDIR}/groovy ${RESOURCES_GOBYWEB_SERVER_SIDE_QUEUE_WRITER} ${QUEUE_WRITER_POSTFIX} "$@"

else
  echo "Status update disabled: $*"
fi
