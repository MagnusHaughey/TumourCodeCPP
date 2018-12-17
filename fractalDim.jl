# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end





#================== Import Python/NumPy modules ===================#
println("")
print("Importing Python/NumPy modules...")
using PyCall
@pyimport numpy as np

print("\rImporting Python/NumPy modules... done.")




#================== Read in import data file from command line ===================#
println("")
print("Reading tumour data...")
path = ARGS[1]
infile = string(path , "/" , ARGS[2] , "_boundaries.csv")
data = np.loadtxt( infile , unpack=true , delimiter="," , skiprows=1 , usecols=(0,1))

# Create x, y, z, GA vectors
data_set = length(data)/2
x1 = fill(0.0 , float_to_int(data_set))
x2 = fill(0.0 , float_to_int(data_set))

# Unpack into separate x, y, z, GA vectors
for i in 1:length(data)
	if ( ((i-1)%2) == 0 ) x1[float_to_int(((i-1)/2)+1)] = data[i] 
	elseif ( ((i-2)%2) == 0 ) x2[float_to_int(((i-2)/2)+1)] = data[i] 
	end
end

print("\rReading tumour data... done.")



#================== Separate different sub-clones ===================#
println("")
print("Separating sub-clones...")

subclones = fill([] , length(x1))
new = 2

for i in 1:length(x1)

	assigned = false

	print("\rSeparating sub-clones... $(float_to_int(i*100.0/length(x1)))%")
	
	# Check if any neighbours of this cell is contained within an established sub-clone 
	for x in -1:1
		for y in -1:1

			if ((x == 0) && (y == 0)) continue end

			for j in 1:length(x1)
				if ( (x1[i] + x , x2[i] + y) in subclones[j] )
					push!(subclones[j] , (x1[i],x2[i]))
					assigned = true
				end
			end	

			if (assigned) continue end

		end
		if (assigned) continue end
	end

	if !(assigned)
		push!(subclones[new] , (x1[i],x2[i]))
		global new += 1
	end

end

print("\rSeparating sub-clones... done.")



#================== Write sub-clone boundary data ===================#
println("")
print("Writing data...")

# Write 2D sliced data
outfile = string(path , "/" , ARGS[2] , "_sep_boundaries.csv" )
open(outfile, "w") do f
	label = 1
	for subclone in subclones
		for point in subclone
			write( f , "$(point[1]),$(point[2]),$label\n")
		end
		label += 1
	end
end

print("\rWriting data... done.")
println("")










