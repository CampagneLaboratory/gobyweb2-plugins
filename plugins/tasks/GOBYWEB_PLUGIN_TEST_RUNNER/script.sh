
export M2_HOME="${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}"
export MAVEN_OPTS="-Xms256m -Xmx512m"

function plugin_task {
    ${FILESET_COMMAND} --has-fileset TEST_CLASSES
    dieUponError "Input test classes are not available"

    TEST_CLASSES_JAR=`${FILESET_COMMAND} --fetch TEST_CLASSES`
    dieUponError "Failed to fetch test classes files ${TEST_CLASSES_JAR}"
    echo ${TEST_CLASSES_JAR}

    mkdir ./source
    cp "${TEST_CLASSES_JAR}" ./source/JarWithTests.jar
    rm -rf ./additionalTests/
    ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml clean
    ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml process-test-resources
    ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn -f ${JOB_DIR}/pom.xml surefire:test -Dtest=${PLUGINS_TASK_GOBYWEB_PLUGIN_TEST_RUNNER_TEST_NAMES}
    # Surefire report directory will be: ./target/plugins-reports
}
