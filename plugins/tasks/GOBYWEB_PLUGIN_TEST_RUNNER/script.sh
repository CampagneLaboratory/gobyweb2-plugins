


function plugin_task {
    mvn clean
    mvn test
    # Surefire report directory will be: ./target/plugins-reports
    mvn surefire:test -Dtest=*

}

plugin_task