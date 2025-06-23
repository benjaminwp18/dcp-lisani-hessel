IFS=$'\n'
search_dir="./input"
for entry in $(ls $search_dir)
do
    ./src/dehazeDCP input/"$entry" results/"$entry" -s 7 -p 0.1 -w 0.7 -r 30 -e 0.0001 -t 0.1 -u -x
    echo $entry
done