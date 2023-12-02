#!/bin/bash
# the loop will read each text file from the current directory
for filename in `ls *.txt`
do
  # Print the text filename before conversion
  echo "Filename before conversion : $filename"
  # Change the extension of the file txt to docx
  mv -- "$filename" "$(basename -- "$filename" .txt).docx"
done
