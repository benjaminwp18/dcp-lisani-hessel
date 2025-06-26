IFS=$'\n'
search_dir="./imgs/input/tclahe"
for entry in $(ls $search_dir)
do
    ./src/dehazeDCP imgs/input/tclahe/"$entry" imgs/output/tclahe_udcp/"$entry" -s 7 -p 0.1 -w 0.95 -r 30 -e 0.0001 -t 0.1 -u
    echo $entry
done