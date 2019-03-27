# 2019-03-27 11:34:17

Quick script to parse article type and year from a set of NXML files. 

```bash
# start by grabbing the names of all directories to a file
ls -1 | grep ajhg > dirs.txt
# remove any old lists
rm list
# for each directory d in the list of directories
for d in `cat dirs.txt`; do
    # write down to a list of files called list any file names that have Type article inside.
    grep -l "TYPE article" $d/*/*nxml >> list;
done
# add headers to the output CSV file
echo "file,article type,article year" > ~/tmp/julia.csv;
# for each file in the list we created above
for f in `cat list`; do
    # store in variable atype the article-type tag and its contents
    # then use sed to remove the tag
    atype=`grep -o 'article-type="[^"]*"' $f | sed "s/article-type=//"`;
    # do the same for the tag year, removing the tag aftewards using sed
    year=`grep -o '<year>[0-9]*</year>' $f | sed 's/<[\/]*year>//g'`;
    # spit out the file name, atype and year variables to the output file
    echo $f,$atype,$year >> ~/tmp/julia.csv;
done
```

Note that we didn't necessarily have to do the first step per directory.
However, when I tried using */*/*nxml, the argument list was too long, so I had
to break it down.