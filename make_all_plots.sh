

allFiles=$1'/*.csv'
allImages=$1'/%03d.csv.png'
outPath=$1'/out.mp4'

# Plot all spatial data
for file in $allFiles
do

	echo $file
	python3 plot_spatial_data.py $file

done


# Combine all images into movie 
ffmpeg -framerate 10 -i $allImages -c:v libx264 -r 30 -pix_fmt yuv420p $outPath
