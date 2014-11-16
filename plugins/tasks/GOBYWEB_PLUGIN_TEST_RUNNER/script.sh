


function plugin_task {
    ${FILESET_COMMAND} --has-fileset TEST_CLASSES
    dieUponError "Input test classes are not available"

    FILES_WITH_LIST=`${FILESET_COMMAND} --fetch TEST_CLASSES`
    dieUponError "Failed to fetch test classes files ${TEST_CLASSES}"
    echo ${TEST_CLASSES}

    mkdir ./source
    cp $FILES_WITH_LIST ./source/JarWithTests.jar
    rm -rf ./additionalTests/
    mvn clean
    mvn test
    # Surefire report directory will be: ./target/plugins-reports
}

plugin_task