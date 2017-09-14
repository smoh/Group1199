#!/bin/bash
set -e

target_dir="submission"
if [ -d "$target_dir" ]; then
  echo "${target_dir} exists; do nothing"
  exit 1;
fi
mkdir $target_dir

list_of_figures=$(grep includegraphics ms.tex | sed 's/.*{\(.*\)}/figures\/\1/')
echo "grabbing following figures\n";
echo $list_of_figures;
for figfile in $list_of_figures; do
  cp $figfile $target_dir;
done

paperfiles=$(git ls-tree --name-only HEAD | grep -v ".gitignore\|Makefile\|figures\|.*.sh")
for f in $paperfiles; do
  echo "copying $f to $target_dir";
  cp $f $target_dir;
done

echo "compress all files in $submission to $submission.tar.gz"
tar -czvf "$target_dir.tar.gz" $target_dir/*;
