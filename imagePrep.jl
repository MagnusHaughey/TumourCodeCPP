


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
infile = ARGS[1]
infile = string("./DATA/" , infile)
data = np.loadtxt( infile , unpack=true , delimiter="," , skiprows=1 , usecols=(0,1,2,4))

# Create x, y, z, GA vectors
data_set = length(data)/4
#println("$data  $data_set")
x = fill(0.0 , float_to_int(data_set))
y = fill(0.0 , float_to_int(data_set))
z = fill(0.0 , float_to_int(data_set))
res = fill(0.0 , float_to_int(data_set))

# Unpack into separate x, y, z, GA vectors
for i in 1:length(data)
	if ( ((i-1)%4) == 0 ) x[float_to_int(((i-1)/4)+1)] = data[i] 
	elseif ( ((i-2)%4) == 0 ) y[float_to_int(((i-2)/4)+1)] = data[i] 
	elseif ( ((i-3)%4) == 0 ) z[float_to_int(((i-3)/4)+1)] = data[i] 
	elseif ( ((i-4)%4) == 0 ) res[float_to_int(((i-4)/4)+1)] = data[i] 
	end
end

print("\rReading tumour data... done.")

#println(x)
#println(y)
#println(z)
#println(res)


#================== Slice 3D data at specified coordinate ===================#
println("")
print("Slicing data...")
slice = ARGS[2]
index = findfirst(isequal('=') , slice)
slice_fraction = parse(Float64 , slice[index+1:end])


if ((slice_fraction < 0.0) || (slice_fraction > 1.0))
	println("Error with 2nd argument. Enter value between 0 and 1. Exiting...")
	exit(0)
end


if (slice[1] == 'x')
	slice_val = findmin(x)[1] + float_to_int( (findmax(x)[1] - findmin(x)[1])*slice_fraction )
elseif (slice[1] == 'y')
	slice_val = findmin(y)[1] + float_to_int( (findmax(y)[1] - findmin(y)[1])*slice_fraction )
elseif (slice[1] == 'z')
	slice_val = findmin(z)[1] + float_to_int( (findmax(z)[1] - findmin(z)[1])*slice_fraction )
else
	println("Error with 2nd argument. Co-ordinate must x, y or z. Exiting...")
	exit(0)
end

slice1 = []
slice2 = []
GAdata = []

if (slice[1] == 'x')
	for i in 1:length(x)
		if (x[i] == slice_val)
			push!(slice1 , y[i])
			push!(slice2 , z[i])
			push!(GAdata , res[i])
		end
	end

elseif (slice[1] == 'y')
	for i in 1:length(y)
		if (y[i] == slice_val)
			push!(slice1 , x[i])
			push!(slice2 , z[i])
			push!(GAdata , res[i])
		end
	end

else
	for i in 1:length(z)
		if (z[i] == slice_val)
			push!(slice1 , x[i])
			push!(slice2 , y[i])
			push!(GAdata , res[i])
		end
	end
end

print("\rSlicing data... done.")


#println(slice1)
#println(slice2)
#println(GAdata)
#println(slice_val)



#================== Algorithm for finding sub-clone boundaries ===================#
println("")
print("Finding boundaries...")

dim = float_to_int(length(GAdata)^(0.5))
boundary = []

for i in 1:length(GAdata)
	if (GAdata[i] != 0)		# find mutated cell in vector

		if ( ((i == 1) || (i == length(GAdata))) && (GAdata[i] == 1.0) )		# If cell has mutation and is on the edge of the tumour
			push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			continue

		elseif ( (i > 1) && (i < length(GAdata)) && (slice1[i-1] != slice1[i+1]) && (GAdata[i] == 1.0) )		# If cell has mutation and is on the edge of the tumour
			push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			continue

		elseif ( slice1[i] == findmax(slice1)[1] )
			push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			continue
		end

		for j in 1:length(slice1)	# check neighbouring cells
			if ( ((abs(slice1[j] - slice1[i]) <= 1.0) && (abs(slice2[j] - slice2[i]) <= 1.0)) && (GAdata[j] < 1.0) )
			    push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			    break
			end

		end
	end
end

print("\rFinding boundaries... done.")




#================== Write sub-clone boundary data ===================#
println("")
print("Writing data...")

# Write 2D sliced data
outfile = string(infile[1:length(infile)-4] , "_" , slice , "_raw.csv" )
open(outfile, "w") do f
	for i in 1:length(GAdata)
		if ( GAdata[i] != -1.0 ) write( f , "$(slice1[i]),$(slice2[i]),$(GAdata[i])\n") end
	end
end

# Write boundary data
outfile = string(infile[1:length(infile)-4] , "_" , slice , "_boundaries.csv" )
open(outfile, "w") do f
	for point in boundary
		write( f , "$(point[1]),$(point[2]),$(point[3])\n")
	end
end

print("\rWriting data... done.")
println("")











