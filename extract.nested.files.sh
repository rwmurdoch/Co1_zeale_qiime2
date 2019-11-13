
#simple command that will grab all files within subfolders in a directory
#("reads" in the given code) and move them to the indicated directory (-t)
#this is especially useful for instances when the Illumina machine provides read files
#nested in respective sample folders

find nested_reads -type f -print0 | xargs -0 mv -t reads
