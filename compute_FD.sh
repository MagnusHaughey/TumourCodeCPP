

number_mutants_file=$1'/number_mutants.dat'

for file in $1/*.csv
do

	#printf '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n'
	#echo $file


	# Find the number of mutant cells 
	number_mutated_cells=$(cat $file | cut -d"," -f3 | grep "1" | wc -l)
	sweep=${file##*/}
	sweep=${sweep%.csv}

	printf "$sweep $number_mutated_cells \n" >> $number_mutants_file



	############ Extract only coordinates of mutated cells for fractal dimension analysis

	boundaries_file=$file'.boundaries.csv'

	printf "Creating: $boundaries_file \n"

	for line in $(cat $file)
	do

		num_mutations=${line##*,}

		if (( num_mutations == 1 ))
		then
			echo $line >> $boundaries_file
		fi
	done


	############ Compute the fractal dimension
	printf "Calculating fractal dimension \n"
	julia -p 4 fractalDimension.jl -q $boundaries_file



done


############ Collect all fractal dimension values
for file in $1/*.csv.boundaries.csv.fractalDimensions/fractalDimensions.csv
do

	all_dimensions=$1/'all_dimensions.dat'

	#echo $file
	sweep=$(echo $file | cut -d"/" -f7)
	sweep=$(echo $sweep | cut -d"." -f1)
	#echo $temp
	#sweep=${temp##*/}
	printf "$sweep $(cat $file)\n" >> $all_dimensions

done



















