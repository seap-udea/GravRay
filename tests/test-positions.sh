#!/bun/bash
# Example:
#    bash tests/test-positions.sh "07/19/2015 00:00:00.000 UTC"

date=$1
if [ "x$date" = "x" ];then echo "You must provide a valid date.";exit 1;fi

scriptname=$(basename $0)
filename="${scriptname%.*}"

# OBJECTS
objects="SUN MERCURY VENUS EARTH MOON MARS_BARYCENTER JUPITER_BARYCENTER SATURN_BARYCENTER URANUS_BARYCENTER NEPTUNE_BARYCENTER"

# TESTING POSITION
echo "Testing position of $objects..."

{
    echo "DATE: $date"
    for object in $objects
    do
	echo "POSITION OF $object"
	./whereisit.exe "$object" "$date" 2> /dev/null
	echo
    done
}|tee scratch/$filename.log

echo "Position has been calculated and stored in scratch/$filename.log"
