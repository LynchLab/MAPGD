read -p "Are you sure you want to replace $1 with $2? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then 
for dir in . ./data_types ./commands ./bam ./io ./mpi ./sql ./stream_tools/ ./testing ./data_conversion ./raw
do
	`sed -i "s/$1/$2/g" $dir/*.h`
	`sed -i "s/$1/$2/g" $dir/*.cc`
	echo replacing $1 with $2 in $dir/*.h
done
fi
