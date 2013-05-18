# This script connects the artifact manager to the GobyWeb environment. It will be source
# before any artifact installation function is executed.
# JOB_DIR is defined by the caller.
if [ -e ${JOB_DIR}/constants.sh ]; then
    . ${JOB_DIR}/constants.sh
else
  touch ${JOB_DIR}/constants.sh
fi

if [ -e ${JOB_DIR}/auto-options.sh ]; then
    . ${JOB_DIR}/auto-options.sh
else
  touch ${JOB_DIR}/auto-options.sh
fi

# environment file to test genome download and index creation.
export SGE_O_WORKDIR=${JOB_DIR}

# This function must exist in an environment script:

function plugin_install_artifact {
    return 0
}