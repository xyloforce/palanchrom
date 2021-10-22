# sort bed
# intersect beds
for file in $(ls $1/*.bed)
do
export filename=$(basename $file)
export filename="${filename%.*}"
if [[ -z "$(ls $2/*sorted.bed)" ]] ; then
    sort -k1,1 -k2,2n $file > $2/$filename".sorted.bed"
else
    echo "already sorted"
fi
done

# intersect : first get source file then intersect each with this one

source=( $(ls $2/*sorted.bed) )
source=${source[0]}
cat $source > $2"/intersected.bed"

for file in $(ls $2/*sorted.bed)
do
export tempFile=$(date +"%N")"_intersectTemp.bed"
bedtools intersect -sorted -a $2"/intersected.bed" -b $file > $tempFile
sort -k1,1 -k2,2n $tempFile > $2"/intersected.bed"
rm $tempFile
done
