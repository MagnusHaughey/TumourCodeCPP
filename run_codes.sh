

total=0
seed=$1
dseed=$2
s=$3
sweeps=$4

successful_simulations='./2D_DATA/maxSize=100000_sweeps='"$sweeps"'_ut=2_BDratio=0_s='"$s"'/successful_simulations.dat'

while :
do
	
	filePath='./2D_DATA/maxSize=100000_sweeps='"$sweeps"'_ut=2_BDratio=0_s='"$s"'/'"$seed"


	# Simulate data
	printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	printf "Simulating: ./2Dmain_homeostasis_withMutation -q $seed $s $sweeps \n"
	./2Dmain_homeostasis_withMutation -q $seed $s $sweeps

	
	# Before analysing, determine if mutated clone expanded -> simulation was 'successful'
	padded_sweeps=$(printf "%03d\n" $sweeps)
	final_number_mutated_cells=$(cat $filePath'/animation/'$padded_sweeps'.csv' | cut -d"," -f3 | grep "1" | wc -l)

	if (( final_number_mutated_cells > 50))
	then
		echo "seed=$seed 1" >> $successful_simulations
	else
		echo "seed=$seed 0" >> $successful_simulations
		seed=$(( seed + dseed ))
		rm -r $filePath
		continue
	fi


	# Compute fractal dimensions
	bash compute_FD.sh $filePath'/animation/'

	seed=$(( seed + dseed ))
	total=$(( total + 1 ))

	# Compress simulated data
	tar czf $filePath'/animation.tar.gz' $filePath'/animation/'
	rm -r $filePath'/animation/'	

	if (( total == 10 ))
	then

		break

	fi

done
