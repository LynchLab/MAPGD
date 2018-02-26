dir=`echo $1 | cut -d '/' -f 1`
a=`echo $1 | cut -d '/' -f 2`
b=`echo $2 | cut -d '/' -f 2`
read -p "Are you sure you want to move and replace $a with $b in $dir? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then 
mv ${1}.h ${2}.h
mv ${1}.cc ${2}.cc
for dir in . ./data_types ./commands ./io ./mpi ./sql ./testing ./data_conversion ./raw
do
	`sed -i "s/$a/$b/g" $dir/*.h`
	`sed -i "s/$a/$b/g" $dir/*.cc`
	echo replacing $a with $b in $dir/*.h
done
fi
