#!/bin/bash

NARGS="$#"
echo "INFO: NARGS= $NARGS"

if [ "$NARGS" -lt 2 ]; then
	echo "ERROR: Invalid number of arguments...see script usage!"
  echo ""
	echo "**************************"
  echo "***     USAGE          ***"
	echo "**************************"
 	echo "$0 [ARGS]"
	echo ""
	echo "=========================="
	echo "==    ARGUMENT LIST     =="
	echo "=========================="
	echo "*** MANDATORY ARGS ***"
	echo "--fileext=[FILE_EXT] - File extension (e.g. txt, root, fits) to be placed in list"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--rootdir=[ROOT_DIR] - Directory where to start searching file to be placed in list [default=pwd]"
	echo "--fileprefix=[FILE_PREFIX] - File prefix filter [default: none]"
	echo "--stripext - Strip file extension before writing to list? [default=no]"
	echo "--strippath - Strip file path before writing to list? [default=no]"
	echo "--searchdir - Search for directories (instead of files) with given prefix and extension [default=no]"
	echo "--recursive - Search recursively down from ROOT_DIR [default=no]"
	echo "--output=[OUTPUTFILE] - Output file name with file list [default=filelist.txt]"
	echo "=========================="
	exit 1
fi

#######################################
##         PARSE ARGS
#######################################
CURRENTDIR=$PWD
FILE_EXT=""
FILE_EXT_GIVEN=false
FILE_PREFIX=""
STRIP_EXT=false
STRIP_PATH=false
ROOT_DIR="$PWD"
OUTPUTFILE="filelist.txt"
SEARCH_TYPE_FLAG="-type f" # search for files
RECURSIVE_FLAG="-maxdepth 1"

for item in $*
do
	case $item in 
		## MANDATORY ##	
		--fileext=*)
    	FILE_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$FILE_EXT" != "" ]; then
				FILE_EXT_GIVEN=true
			fi
    ;;
		## OPTIONAL ##	
		--rootdir=*)
    	ROOT_DIR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;	
		--fileprefix=*)
    	FILE_PREFIX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--output=*)
    	OUTPUTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		--stripext*)
    	STRIP_EXT=true
    ;;
		--strippath*)
    	STRIP_PATH=true
    ;;
		--recursive*)
			RECURSIVE_FLAG=" "
    ;;
		--searchdir*)
    	SEARCH_TYPE_FLAG="-type d" # search for directories
    ;;

    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done



#######################################
##         PROCESS FILENAME
#######################################
process_filename(){
	local filename=$1
	local strippath=$2
	local stripext=$3
	local fileext=$4
	local outputfile=$5
	
	# Strip path?
	if [ $strippath = true ]; then
		filename_base=$(basename $filename)
		filename=$filename_base
	fi
		
	# Strip extension?
	if [ $stripext = true ]; then
		filename_noext="${filename//.$fileext/}"
		filename=$filename_noext
	fi

	# Print processed filename to file list	
	echo $filename
  echo $filename >> $outputfile
}



## Process file/directory found by find command
FIND_CMD="find $ROOT_DIR $RECURSIVE_FLAG $SEARCH_TYPE_FLAG -name "'"'"$FILE_PREFIX*.$FILE_EXT"'"'" | sort --version-sort"
echo "Executing find command: $FIND_CMD"

#find $ROOT_DIR "$RECURSIVE_FLAG" "$SEARCH_TYPE_FLAG" -name "$FILE_PREFIX*.$FILE_EXT"  | while read item; do
eval $FIND_CMD | while read item; do
	process_filename "$item" "$STRIP_PATH" "$STRIP_EXT" "$FILE_EXT" "$CURRENTDIR/$OUTPUTFILE"
done



