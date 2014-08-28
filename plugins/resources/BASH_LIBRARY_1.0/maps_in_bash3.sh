####################
# This set of functions simulated the behavior of associative arrays (or maps) in Bash 3.x.
# Associative arrays are available only on Bash 4+
#
# Sample usage:
#
# to store some elements in the map:
#
# put "my-map-name" "key1" "value1"
# put "my-map-name" "key2" "value2"
# put "my-map-name" "key3" "value3"
#
# to retrieve the value associated to a key in the map:
#
# get "my-map-name" "key2"
# echo $value (will print value2)
#
# to print all the elements in the map:
#
# getKeySet "my-map-name"
# for KEY in $keySet
# do
#  get "my-map-name" $KEY
#  echo $KEY -> $value
# done
#
# the loop will print:
# key1 -> value1
# key2 -> value2
# key3 -> value3
####################

# populate a map starting from the keys and values in the given file
# the format of each line in the file must be Key=Value
populateMapFromPropertiesFile() {
    mapName=$1; file=$2;

    while read line; do
        line=`echo $line | awk '{split($0,a, /[#!]/); print a[1]}'`
        read -rd '' line <<< "$line" # trim decommented line
        read -rd '' line <<< "$line" # trim decommented line
       if [ -n "$line" ]; then #skip empty
            #echo "[$line]"
            splitidx=`expr index "$line" '[=:]'` # get the first index of the equals sign or colon
            if [ $splitidx = 0 ]; then
                key="$line"
                val=""
            else
                key=${line:0:$splitidx-1} # string before assignment operator
                val=${line:$splitidx} # string after assignment operator
            fi
        fi
        read -rd '' key <<< "$key" # trim string
        key=${key//./|} # nasty ugly hack for associative arrays not being allowed to have periods in the key
        read -rd '' val <<< "$val" #trim string
        put $mapName $key $val
    done < $file
}

# store a new key in a map, if the map is simulated, it doen not need to be created
# parameters: mapname, key, value
put() {
    if [ "$#" != 3 ]; then exit 1; fi
    mapName=$1; key=$2; value=`echo $3 | sed -e "s/ /:SP:/g"`
    eval map="\"\$$mapName\""
    map="`echo "$map" | sed -e "s/--$key=[^ ]*//g"` --$key=$value"
    eval $mapName="\"$map\""
}

# get the value associated to the given key in the map
# parameters: mapname, key
# return: the value of the key in the $value variable
get() {
    mapName=$1; key=$2; valueFound="false"

    eval map=\$$mapName

    for keyValuePair in ${map};
    do
        case "$keyValuePair" in
            --$key=*) value=`echo "$keyValuePair" | sed -e 's/^[^=]*=//'`
                      valueFound="true"
        esac
        if [ "$valueFound" == "true" ]; then break; fi
    done
    value=`echo $value | sed -e "s/:SP:/ /g"`
}

# get all the keys in the map
# parameters: mapname
# return: the list of keys is returned in the $keySet variable
getKeySet() {
    if [ "$#" != 1 ];
    then
        exit 1;
    fi

    mapName=$1;

    eval map="\"\$$mapName\""

    keySet=`
           echo $map |
           sed -e "s/=[^ ]*//g" -e "s/\([ ]*\)--/\1/g"
          `
}
