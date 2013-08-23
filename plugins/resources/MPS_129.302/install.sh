# Installation script for MPS 3.0
MAVEN_CENTRAL_REPO_URL=http://repo1.maven.org/maven2/
CAMPAGNELAB_RELEASE_REPO_URL=http://repository.campagnelab.org/artifactory/CampagneLab/
CAMPAGNELAB_SNAPSHOT_REPO_URL=http://repository.campagnelab.org/artifactory/CampagneLab-SNAPSHOT/
REPO_USER=downloader
REPO_PASSWORD=password goes here

function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'DISTRIBUTION' )
            VERSION="3.0"
            BUILD="RC3-129.302"
            ${RESOURCES_FETCH_URL_SCRIPT} http://download.jetbrains.com/mps/30/MPS-${VERSION}-${BUILD}.tar.gz MPS-${VERSION}-${BUILD}.tar.gz
            gzip -c -d MPS-${VERSION}-${BUILD}.tar.gz |tar -xvf -
            cp -r MPS\ ${VERSION}/* ${installation_path}/
            return 0
            ;;
         'SUPPORT_LIBS' )
                download_wildcard
                download_pluginsSDK
                download_stepslogger
                download_filesets
                download_artifacts
                download_ovl
                download_runtime_support
                download_groovy_all
                download_fastutil
                download_scala_library
                download_common_io
                return 0;
            ;;
        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}

function download_wildcard {
   ${RESOURCES_FETCH_URL_SCRIPT} http://wildcard.googlecode.com/files/wildcard-1.03.zip wildcard-1.03.zip
   unzip wildcard-1.03.zip
   cp wildcard-1.03/wildcard-1.03.jar ${installation_path}/
}

function download_pluginsSDK {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_SNAPSHOT_REPO_URL} \
   org/campagnelab/gobyweb/plugins/2.3-SNAPSHOT/plugins-2.3-20130821.191554-33.jar ${installation_path}/plugins-2.3-SNAPSHOT.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_stepslogger {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_SNAPSHOT_REPO_URL} \
   org/campagnelab/org.campagnelab.stepslogger/1.1-SNAPSHOT/org.campagnelab.stepslogger-1.1-20130731.135018-1.jar ${installation_path}/org.campagnelab.stepslogger-1.1-SNAPSHOT.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_filesets {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_RELEASE_REPO_URL} \
   org/campagnelab/gobyweb/filesets/1.1/filesets-1.1.jar ${installation_path}/filesets-1.1.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_artifacts {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_RELEASE_REPO_URL} \
   org/campagnelab/gobyweb/artifacts/2.2.3/artifacts-2.2.3.jar ${installation_path}/artifacts-2.2.3.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_ovl {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_RELEASE_REPO_URL} \
   org/campagnelab/gobyweb/option-validation-language/2.1.0/option-validation-language-2.1.0.jar ${installation_path}/option-validation-language-2.1.0.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_runtime_support {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${CAMPAGNELAB_SNAPSHOT_REPO_URL} \
   org/campagnelab/nyosh/nyosh-runtime-support/1.0.3-SNAPSHOT/nyosh-runtime-support-1.0.3-20130731.211633-1.jar ${installation_path}/nyosh-runtime-support-1.0.3-SNAPSHOT.jar \
   ${REPO_USER} ${REPO_PASSWORD}
}

function download_groovy_all {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${MAVEN_CENTRAL_REPO_URL} \
   org/codehaus/groovy/groovy-all/1.8.6/groovy-all-1.8.6.jar ${installation_path}/groovy-all-1.8.6.jar \
   '' ''
}

function download_fastutil {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${MAVEN_CENTRAL_REPO_URL} \
   it/unimi/dsi/fastutil/6.4.4/fastutil-6.4.4.jar ${installation_path}/fastutil-6.4.4.jar \
   '' ''
}

function download_scala_library {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${MAVEN_CENTRAL_REPO_URL} \
   org/scala-lang/scala-library/2.9.2/scala-library-2.9.2.jar ${installation_path}/scala-library-2.9.2.jar \
   '' ''
}

function download_common_io {
   ${RESOURCES_MAVEN_ARTIFACTS_DOWNLOADER_RUN_DOWNLOADER} ${MAVEN_CENTRAL_REPO_URL} \
   commons-io/commons-io/2.4/commons-io-2.4.jar ${installation_path}/commons-io-2.4.jar \
   '' ''
}