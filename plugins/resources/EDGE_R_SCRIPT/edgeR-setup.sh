#!/bin/bash
# Execute the samtools executable from the artifact repository:
. ${JOB_DIR}/artifacts.sh
expose_artifact_environment_variables

${RESOURCES_ARTIFACTS_EDGER_BINARIES}/setup.sh
