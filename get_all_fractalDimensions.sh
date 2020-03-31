


tar_archives=$1'/*/animation.tar.gz'
num_sweeps=$2
tar_ext=".tar.gz"

for archive in $tar_archives
do

	# Remove any double slashes from file path. This seems to cause an error when reading data file in tar archive
	#archive=$(readlink -m $archive)
	archive=$(echo $archive | sed s#//*#/#g)
	archive_noExtension=${archive%"$tar_ext"}

	for sweep in $(seq 0 $num_sweeps)
	#for file in $(tar -ztvf $archive | rev | cut -d" " -f1 | rev | grep "fractalDimensions.csv")
	do

		padded_sweep=$(printf "%03d\n" $sweep)
		fractalData=$archive_noExtension'/'$padded_sweep'.csv.boundaries.csv.fractalDimensions/fractalDimensions.csv'
		#echo $sweep

		
		#echo $archive
		#echo $fractalData
		#printf "tar -axf $archive $fractalData -O"
		fractalDimensionValue=$(tar -axf $archive $fractalData -O)

		printf "$sweep $fractalDimensionValue \n"

	done

done
