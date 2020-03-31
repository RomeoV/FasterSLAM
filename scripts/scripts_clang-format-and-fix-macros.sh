#!/bin/bash

# If the first argument is just '-', do not replace the files
# instead print the result to stdout.

sed_i="-i .bak"
testername=/tmp/sed_tester.$$.tmp
touch $testername
rm -f $testername.bak
sed ${sed_i} -e '' $testername > /dev/null 2>&1
result=$?
if [ -e $testername.bak ] ; then
  sed_i='-i ""'
else
  sed_i='-i '
fi

macros=(
	__host__
	__device__
	VECCORE_ATT_HOST
	VECCORE_ATT_DEVICE
	VECCORE_ATT_HOST_DEVICE
	VECCORE_FORCE_INLINE
	VECCORE_FORCE_NOINLINE
	VECGEOM_CUDA_HEADER_HOST
	VECGEOM_CUDA_HEADER_DEVICE
	VECGEOM_CUDA_HEADER_BOTH
	VECGEOM_FORCE_INLINE
)

macroShorts=(
# Obsolete
    "VECGEOM_FORCE_INLINE[[:space:]]*VECGEOM_CUDA_HEADER_BOTH @GIGB \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECGEOM_CUDA_HEADER_BOTH[[:space:]]*VECGEOM_FORCE_INLINE @GIGB \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECGEOM_CUDA_HEADER_BOTH @GB \1VECCORE_ATT_HOST_DEVICE\1"
    "VECGEOM_FORCE_INLINE[[:space:]]*VECGEOM_CUDA_HEADER_DEVICE @GIGD \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECGEOM_CUDA_HEADER_DEVICE[[:space:]]*VECGEOM_FORCE_INLINE @GIGD \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECGEOM_CUDA_HEADER_DEVICE @GD \1VECCORE_ATT_DEVICE\1"
    "VECGEOM_FORCE_INLINE[[:space:]]*VECGEOM_CUDA_HEADER_HOST @GIGH \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST\1"
    "VECGEOM_CUDA_HEADER_HOST[[:space:]]*VECGEOM_FORCE_INLINE @GIGH \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST\1"
    "VECGEOM_CUDA_HEADER_HOST @GH \1VECCORE_ATT_HOST\1"

# Mixed.
    "VECGEOM_FORCE_INLINE[[:space:]]*VECCORE_ATT_HOST_DEVICE @GICB \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECCORE_ATT_HOST_DEVICE[[:space:]]*VECGEOM_FORCE_INLINE @GICB \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECGEOM_FORCE_INLINE[[:space:]]*VECCORE_ATT_DEVICE @GICD \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECCORE_ATT_DEVICE[[:space:]]*VECGEOM_FORCE_INLINE @GICD \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECGEOM_FORCE_INLINE[[:space:]]*VECCORE_ATT_HOST @GICH \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST\1"
    "VECCORE_ATT_HOST[[:space:]]*VECGEOM_FORCE_INLINE @GICH \1VECGEOM_FORCE_INLINE\1VECCORE_ATT_HOST\1"

    "VECGEOM_FORCE_INLINE @GI \1VECGEOM_FORCE_INLINE\1"

# VecCore only
    "VECCORE_FORCE_INLINE[[:space:]]*VECCORE_ATT_HOST_DEVICE @CICB \1VECCORE_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECCORE_ATT_HOST_DEVICE[[:space:]]*VECCORE_FORCE_INLINE @CICB \1VECCORE_FORCE_INLINE\1VECCORE_ATT_HOST_DEVICE\1"
    "VECCORE_ATT_HOST_DEVICE @CB \1VECCORE_ATT_HOST_DEVICE\1"
    "VECCORE_FORCE_INLINE[[:space:]]*VECCORE_ATT_DEVICE @CICD \1VECCORE_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECCORE_ATT_DEVICE[[:space:]]*VECCORE_FORCE_INLINE @CICD \1VECCORE_FORCE_INLINE\1VECCORE_ATT_DEVICE\1"
    "VECCORE_ATT_DEVICE @CD \1VECCORE_ATT_DEVICE\1"
    "VECCORE_FORCE_INLINE[[:space:]]*VECCORE_ATT_HOST @CICH \1VECCORE_FORCE_INLINE\1VECCORE_ATT_HOST\1"
    "VECCORE_ATT_HOST[[:space:]]*VECCORE_FORCE_INLINE @CICH \1VECCORE_FORCE_INLINE\1VECCORE_ATT_HOST\1"
    "VECCORE_ATT_HOST @CH \1VECCORE_ATT_HOST\1"
    "VECCORE_FORCE_INLINE @CI \1VECCORE_FORCE_INLINE\1"
    "VECCORE_FORCE_NOINLINE @CNI \1VECCORE_FORCE_NOINLINE\1"
    "__host__ @H \1__host__\1"
    "__device__ @D \1__device__\1"
    "__host__ __device__ @B \1__host__\1__device__\1"
    "__device__ __host__ @B \1__host__\1__device__\1"
)

patternFrom=""
patternTo=""

for macroInfo in "${macroShorts[@]}" ; do
    macro=${macroInfo%% *}
    values=${macroInfo#* }
    short=${values%% *}
    newPattern=${values#* }
    # printf "%s switch to %s then %s\n" "$macro" "$short" "$newPattern"

    patternFrom="${patternFrom}
s/\(\n[[:blank:]]*template[^;{(]*\)>[[:space:]]*${macro}/\1${short}>/g"

    patternTo="${patternTo}
s/ *${short}>\([[:space:]\n]*\)/>${newPattern}/g"
done

if [ "-" = "$1" ] ; then
  shift;

  for file in $@; do
    cat $file | eval sed -n -e '"
      # if the first line copy the pattern to the hold buffer
      1h
      # if not the first line then append the pattern to the hold buffer
      1!H
      # if the last line then ...
      $ {
        # copy from the hold to the pattern buffer
        g
        # do the search and replace
        ${patternFrom}
        # print
        p
      }"' | clang-format  | clang-format | eval sed -n -e '"
      # if the first line copy the pattern to the hold buffer
      1h
      # if not the first line then append the pattern to the hold buffer
      1!H
      # if the last line then ...
      $ {
        # copy from the hold to the pattern buffer
        g
        # do the search and replace
        ${patternTo}
        # print
        p
      }"'
  done
else
# Shorten the macro and move them out of the way so that they have
# no effect on line length calculation.
  for file in $@; do
    eval sed  ${sed_i} -n -e '"
      # if the first line copy the pattern to the hold buffer
      1h
      # if not the first line then append the pattern to the hold buffer
      1!H
      # if the last line then ...
      $ {
        # copy from the hold to the pattern buffer
        g
        # do the search and replace
        ${patternFrom}
        # print
        p
      }"' '"${file}"'
  done

  clang-format -i "$@"
  # Run clang-format a 2nd time, this stabilizes some of the comment positioning.
  clang-format -i "$@"

  # Put back the macros.
  for file in $@; do
    eval sed   ${sed_i} -n -e '"
      # if the first line copy the pattern to the hold buffer
      1h
      # if not the first line then append the pattern to the hold buffer
      1!H
      # if the last line then ...
      $ {
        # copy from the hold to the pattern buffer
        g
        # do the search and replace
        ${patternTo}
        # print
        p
      }"' '"${file}"'
  done
fi