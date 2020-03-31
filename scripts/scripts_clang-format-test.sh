#!/bin/bash

#
# This is running a clang-format test
# by doing a filtering step and then analysing
# the result of applying ./scripts/clang-format-and-fix-macros.sh
#

./scripts/clang-format-apply.sh
res=$?

if [ $res -eq 0 ] ; then 
   # cleanup changes in git
   git reset HEAD --hard
fi

exit $res
