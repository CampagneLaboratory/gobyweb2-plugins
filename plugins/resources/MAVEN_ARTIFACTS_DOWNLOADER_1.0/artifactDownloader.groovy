/**
 * This groovy script downloads artifacts from the CampagneLab's Maven Repository using the REST API exposed by artifactory
 * @author manuele
 */
@Grab(group='org.codehaus.groovy.modules.http-builder', module='http-builder', version='0.6' )
import groovyx.net.http.RESTClient
@Grab(group='org.apache.httpcomponents', module='httpcore', version='4.3')
import org.apache.http.entity.FileEntity

//input script parameters
String REPO = args[0]  //the MAVEN repository URL (e.g. http://repository.campagnelab.org/artifactory/CampagneLab-SNAPSHOT/)
String ARTIFACT_PATH = args[1] //the path of the artifact as stored in the REPO (e.g. org/campagnelab/gobyweb/plugins/2.3-SNAPSHOT/plugins-2.3-20130821.191554-33-sdk.jar)
String LOCAL_PATH = args[2] //the path of the file where to store the artifact  (e.g. plugins-2.3-SNAPSHOT.jar)
String USER = args[3] //user of REPO with read permissions
String PASSWORD = args[4]

/**
 * print the byte stream to the output file
 * @param data
 */
def toFile(ByteArrayInputStream data, String file) {

    FileOutputStream fos = new FileOutputStream(file)
    int size = data.available()
    byte[] bytes = new byte[size]
    data.read(bytes, 0, size)
    for (int i = 0; i < size;)
        fos.write(bytes[i++]&0xff)
    fos.close();
}

def restClient = new RESTClient(REPO)
restClient.auth.basic(USER, PASSWORD)
def response = restClient.get(path: ARTIFACT_PATH, requestContentType: 'application/java-archive')
assert response.status == 200
toFile(response.getData(), LOCAL_PATH)

