#!/bin/bash

export TARGET="/media/tom/cfd04a45-d8e3-4944-a080-6cfef1543df8/home/tom/Documents/writing/lmdPaper/data"
export DESTINATION=data

# make directories
find $TARGET -type d -print0 | xargs -0 bash -c 'for DIR in "$@"; 
do
  mkdir -p $DESTINATION${DIR#$TARGET}
  done' -

# make links
find $TARGET -type f -name "*.Rds" -print0 | xargs -0 bash -c 'for file in "$@"; 
do
  ln -s "$file" "$DESTINATION"${file#$TARGET}
  done' -