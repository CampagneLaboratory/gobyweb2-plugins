
    export M2_HOME="${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}"
    export MAVEN_OPTS="-Xms512m -Xmx1g"

    function plugin_task {
        ${FILESET_COMMAND} --has-fileset TEST_CLASSES
        dieUponError "Input test classes are not available"

        TEST_CLASSES_JAR=`${FILESET_COMMAND} --fetch TEST_CLASSES`
        dieUponError "Failed to fetch test classes files ${TEST_CLASSES_JAR}"
        echo ${TEST_CLASSES_JAR}
        mkdir ${JOB_DIR}/source
        cp "${TEST_CLASSES_JAR}" ${JOB_DIR}/source/JarWithTests.jar
        rm -rf ./additionalTests/
        ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml clean
        ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml process-test-resources
        ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml surefire:test -Dtest=${PLUGINS_TASK_GOBYWEB_PLUGIN_TEST_RUNNER_TEST_NAMES}
        # Surefire report directory will be: ./target/plugins-reports
        REPORT=`cat ${JOB_DIR}/target/plugins-reports/*.txt`
        info "${REPORT}" 'post_process'
        if [ -n "${PLUGINS_TASK_GOBYWEB_PLUGIN_TEST_RUNNER_COPY_BACK_LOCATION}" ]; then
          scp -r ${JOB_DIR}/target/plugins-reports  ${PLUGINS_TASK_GOBYWEB_PLUGIN_TEST_RUNNER_COPY_BACK_LOCATION}
          info "Test results copied at ${PLUGINS_TASK_GOBYWEB_PLUGIN_TEST_RUNNER_COPY_BACK_LOCATION}" 'post_process'
        fi
    }