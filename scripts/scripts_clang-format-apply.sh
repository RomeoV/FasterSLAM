#!/bin/bash

#
# This is apllying clang-format 
# by doing a filtering step and then 
#  applying ./scripts/clang-format-and-fix-macros.sh
#

# check that we are in a clean state in order to prevent accidential
# changes
if [ ! -z "$(git status --untracked-files=no  --porcelain)" ]; then 
  echo "Script must be applied on a clean git state"
  exit 1
fi


# Retrieve list of files that were changed in source branch
# with respect to master (target branch)
tb="master"
# we disregard the test/static_analysis directory since this contains LLVM code + headers
filelist=`git diff ${tb}... --name-only | grep -v -i -e linkdef | grep -v -e nudy -e cmsToyGV -e "examples/physics/[^/]*/Geant4"`

# function to check if C++ file (based on suffix)
# can probably be done much shorter
function checkCPP(){
    if [[ $1 == *.cc ]];then
	return 0
    elif [[ $1 == *.cpp ]];then
	return 0
    elif [[ $1 == *.cxx ]];then
	return 0
    elif [[ $1 == *.C ]];then
	return 0
    elif [[ $1 == *.c++ ]];then
	return 0
    elif [[ $1 == *.c ]];then
	return 0
    elif [[ $1 == *.CPP ]];then
	return 0
	# header files
    elif [[ $1 == *.h ]];then
	return 0
    elif [[ $1 == *.hpp ]];then
	return 0
    elif [[ $1 == *.hh ]];then
	return 0
    elif [[ $1 == *.icc ]];then
	return 0
    fi
    return 1
}

echo
echo "Checking formatting using the following clang-format version:"
clang-format --version
echo 

# check list of files
for f in $filelist; do
    if checkCPP $f; then
	echo "CHECKING MATCHING FILE ${f}"
	# apply the clang-format script
	./scripts/clang-format-and-fix-macros.sh ${f}
    fi
done

# check if something was modified
notcorrectlist=`git status --porcelain | grep '^ M' | cut -c4-`
# if nothing changed ok
if [[ -z $notcorrectlist ]]; then
  # send a negative message to gitlab
  echo "Excellent. **VERY GOOD FORMATTING!** :thumbsup:"
  exit 0;
else
  echo "The following files have clang-format problems (showing patches)";
  git diff $notcorrectlist
fi

exit 1



