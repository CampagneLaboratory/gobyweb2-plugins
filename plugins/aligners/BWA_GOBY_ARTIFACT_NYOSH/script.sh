# This is the only function that aligners need to implement.
# Parameters:
#   $1: a temporary filename
#   $2: the basename that should be used to store the sorted alignment

function plugin_align {
  #sample parameters reading
  OUTPUT=$1
  BASENAME=$2
  #invoke the model through the script generated by RunMpsScript
  . ${JOB_DIR}/run_model.sh plugin_align ${OUTPUT} ${BASENAME}

}
 