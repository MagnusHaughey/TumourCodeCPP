


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
infile = string(path , "/tumour.csv")
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


xmin = findmin(x)[1]
ymin = findmin(y)[1]
zmin = findmin(z)[1]

xmax = findmax(x)[1]
ymax = findmax(y)[1]
zmax = findmax(z)[1]



if (slice[1] == 'x')
	slice_val = xmin + float_to_int( (xmax - xmin)*slice_fraction )
	sliced_tumour = fill(-1.0 , (float_to_int( ymax - ymin + 1  ),float_to_int( zmax - zmin + 1 )) )

elseif (slice[1] == 'y')
	slice_val = ymin + float_to_int( (ymax - ymin)*slice_fraction )
	sliced_tumour = fill(-1.0 , (float_to_int( xmax - xmin + 1  ),float_to_int( zmax - zmin + 1 )) )

elseif (slice[1] == 'z')
	slice_val = zmin + float_to_int( (zmax - zmin)*slice_fraction )
	sliced_tumour = fill(-1.0 , (float_to_int( xmax - xmin + 1  ),float_to_int( ymax - ymin + 1 )) )

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
			sliced_tumour[float_to_int(y[i] - ymin) + 1 , float_to_int(z[i] - zmin) + 1] = res[i]
			#push!(slice1 , y[i])
			#push!(slice2 , z[i])
			#push!(GAdata , res[i])
			#print("\rSlicing data... $(i*100.0/(dim^2))%")
		end
	end

elseif (slice[1] == 'y')
	for i in 1:length(y)
		if (y[i] == slice_val)
			sliced_tumour[float_to_int(x[i] - xmin) + 1 , float_to_int(z[i] - zmin) + 1] = res[i]
			#push!(slice1 , x[i])
			#push!(slice2 , z[i])
			#push!(GAdata , res[i])
		end
	end

else
	for i in 1:length(z)
		if (z[i] == slice_val)
			sliced_tumour[float_to_int(x[i] - xmin) + 1 , float_to_int(y[i] - ymin) + 1] = res[i]
			#push!(slice1 , x[i])
			#push!(slice2 , y[i])
			#push!(GAdata , res[i])
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


boundary = []


for i in 1:size(sliced_tumour)[1]
	for j in 1:size(sliced_tumour)[2]

		if (sliced_tumour[i,j] > 0.0)

			# Check 8 nearest neighbours
			if ( (i > 1) && (j > 1) && (i < size(sliced_tumour)[1]) && (j < size(sliced_tumour)[2]) )
				for x in -1:1
					for y in -1:1

						if ((x == 0) && (y == 0)) 
							#println("$(i+x) , $(j+y) ***") 
							continue 
						end

						#println("$(i+x) , $(j+y)")
						if ( (sliced_tumour[i+x , j+y] <= 0.0) && !((i,j,sliced_tumour[i,j]) in boundary) )
							push!(boundary , (i , j , sliced_tumour[i,j]))
							continue
						end
					end
				end

			elseif ((i > 1) || (j > 1) || (i < size(sliced_tumour)[1]) || (j < size(sliced_tumour)[2]))
				push!(boundary , (i , j , sliced_tumour[i,j]))
				continue
			end

		end
	end
end




#=
for i in 1:length(GAdata)
	if (GAdata[i] > 0.0)		# find mutated cell in vector
		print("---------------------------------\n")
		println("$(slice1[i]) , $(slice2[i])")

		# Check 8 nearest neighbours
		if ( (i > 1) && (i < length(GAdata)) )
			if ( (GAdata[i-1] <= 0.0) || (GAdata[i+1] <= 0.0) )

				println("$(slice1[i-1]) , $(slice2[i-1])")
				println("$(slice1[i+1]) , $(slice2[i+1])")
				
				push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
				continue	
			
			end
		end

		if ( (i > (dim + 1)) && (i < length(GAdata) - dim) )
			for j in -1:1

				if ((GAdata[i-dim+j] <= 0.0) || (GAdata[i+dim+j] <= 0.0))
					println("$(slice1[i-dim+j]) , $(slice2[i-dim+j])")
					println("$(slice1[i+dim+j]) , $(slice2[i+dim+j])")

					push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
					continue
				
				end
			end
		end


		# Take care of boundary cells
		if ( (i == 1) || (i == length(GAdata)) )

			push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			continue

		elseif !(i > (dim + 1)) || !(i < length(GAdata) - dim)

			push!(boundary , (slice1[i] , slice2[i] , GAdata[i]))
			continue

		end

#=
		if ( ((i == 1) || (i == length(GAdata))) && (GAdata[i] >= 1.0) )		# If cell has mutation and is on the edge of the tumour
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


=#
	end
end
=#
print("\rFinding boundaries... done.")




#================== Write sub-clone boundary data ===================#
println("")
print("Writing data...")

# Write 2D sliced data
outfile = string(path , "/" , slice , "_raw.csv" )
open(outfile, "w") do f
	#for i in 1:length(GAdata)
	#	if ( GAdata[i] != -1.0 ) write( f , "$(slice1[i]),$(slice2[i]),$(GAdata[i])\n") end
	#end
	for i in 1:size(sliced_tumour)[1]
		for j in 1:size(sliced_tumour)[2]
			if ( sliced_tumour[i,j] != -1.0 ) write( f , "$i,$j,$(sliced_tumour[i,j])\n") end
		end
	end
end

# Write boundary data
outfile = string(path , "/" , slice , "_boundaries.csv" )
open(outfile, "w") do f
	for point in boundary
		write( f , "$(point[1]),$(point[2]),$(point[3])\n")
	end
end

print("\rWriting data... done.")
println("")











