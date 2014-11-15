


function plugin_task {
    mvn clean
    mvn test
    # Surefire report directory will be: ./target/plugins-reports
}

plugin_task