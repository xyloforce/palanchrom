touch $(basename $1)".fna"
cumulated=$(basename $1)".fna"
for file in $(ls $1/*.fna)
do
export filename=$(basename $file)
export filename="${filename%.*}"
echo ">"$filename >> $cumulated
sed '1d' $file >> $cumulated
done

bedtools getfasta -fi $cumulated -bed $2 > $3
