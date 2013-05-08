
function calculateMD5 {
 md5url=$1
 RESULT=

 #the following should work on linux
 if [[ -x `which md5sum` ]]; then
     RESULT=`echo ${md5url} | md5sum |cut -d " " -f 1`
 fi

 #the following should work on Mac OS X
 if [[ -x `which md5` ]]; then
     RESULT=`echo ${md5url} | md5`
 fi

}