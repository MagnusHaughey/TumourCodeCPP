


# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end




function sort_labels( subclones )


	for i in 1:findmax(subclones)[1]

		next_label = 0

		if !(i in subclones)

			for q in (i+1):(findmax(subclones)[1])
				if (q in subclones) 
					next_label = q 
					break
				end 
			end

			if (next_label != 0)
				for j in 1:size(subclones)[1]
					for k in 1:size(subclones)[2]
						if (subclones[j,k] == next_label) 
							subclones[j,k] = i 
							#println("Shifted $next_label to $i")
						end
					end
				end
			end
		end
	end

end



#================== Import Python/NumPy modules ===================#
println("")
print("Importing Python/NumPy modules...")
using PyCall
using Pandas
@pyimport numpy as np
@pyimport numpy.core.defchararray as char

print("\rImporting Python/NumPy modules... done.")




#================== Read in import data file from command line ===================#
println("")
print("Reading tumour data...")
path = ARGS[1]
infile = string(path , "/tumour.csv")
#data = np.loadtxt( infile , unpack=true , delimiter="," , skiprows=1 , usecols=(0,1,2,4))


# Determine if data was simulated with or without cell death
death = path[char.find(path , "death")[1] + 7]


# Read simulated tumour data file
data = Array(read_csv(infile, usecols=(0,1,2,4)))



print("\rReading tumour data... done.")



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




# Find max and min x, y and z values in data array
xmin = data[1,1]
xmax = data[1,1]

ymin = data[1,2]
ymax = data[1,2]

zmin = data[1,3]
zmax = data[1,3]


for i in 1:size(data)[1]

	if (data[i,1] < xmin) global xmin = data[i,1] end
	if (data[i,1] > xmax) global xmax = data[i,1] end

	if (data[i,2] < ymin) global ymin = data[i,2] end
	if (data[i,2] > ymax) global ymax = data[i,2] end

	if (data[i,3] < zmin) global zmin = data[i,3] end
	if (data[i,3] > zmax) global zmax = data[i,3] end

end



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
	for i in 1:size(data)[1]
		if (data[i,1] == slice_val)
			sliced_tumour[float_to_int(data[i,2] - ymin) + 1 , float_to_int(data[i,3] - zmin) + 1] = data[i,4]
		end
	end

elseif (slice[1] == 'y')
	for i in 1:size(data)[1]
		if (data[i,2] == slice_val)
			sliced_tumour[float_to_int(data[i,1] - xmin) + 1 , float_to_int(data[i,3] - zmin) + 1] = data[i,4]
		end
	end

else
	for i in 1:size(data)[1]
		if (data[i,3] == slice_val)
			sliced_tumour[float_to_int(data[i,1] - xmin) + 1 , float_to_int(data[i,2] - ymin) + 1] = data[i,4]
		end
	end
end

print("\rSlicing data... done.")




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


print("\rFinding boundaries... done.")



#================== Separate different sub-clones ===================#
println("")
print("Separating sub-clones...")

x1 = []
x2 = []
for point in boundary
	push!(x1, point[1])
	push!(x2, point[2])
end
	
sizeX = float_to_int(findmax(x1)[1] - findmin(x1)[1]) + 1
sizeY = float_to_int(findmax(x2)[1] - findmin(x2)[1]) + 1
subclones = fill(0.0 , (sizeX,sizeY))

#print(sizeX)
#println("\n")

minX = float_to_int(findmin(x1)[1])
minY = float_to_int(findmin(x2)[1])

#println("$minX , $minY")
new = 0

for i in 1:length(x1)

	#println("--------------------------------------")
	#println("$(x1[i]) , $(x2[i])")

	assigned = false
	coordX = float_to_int(x1[i])-minX+1
	coordY = float_to_int(x2[i])-minY+1
	#println(" ($(x1[i]) , $(x2[i])) --- ($(coordX) , $(coordY)) ")
	
	# Check if any neighbours of this cell is contained within an established sub-clone 
	for x in -1:1

		#if (assigned) break end	

		for y in -1:1

			#if (assigned) break end

			if ((x == 0) && (y == 0)) continue end

			# Take care of boundaries
			if ( (coordX+x < 1) || (coordX+x > sizeX) || (coordY+y < 1) || (coordY+y > sizeY) )
				continue
			end

			# Check neighbours
			if ( subclones[coordX+x , coordY+y] != 0.0 )

				if !(assigned)
					subclones[coordX , coordY] = subclones[coordX+x , coordY+y]
					assigned = true
				
				# Join any subclones up if they are really part of the same subclone
				elseif ( (subclones[coordX , coordY] != subclones[coordX+x , coordY+y]) && (subclones[coordX+x , coordY+y] != 0) )
					#println("CLASH BETWEEN ($coordX , $coordY)=$(subclones[coordX , coordY]) and ($(coordX+x) , $(coordY+y))=$(subclones[coordX+x , coordY+y])\n")
					for j in 1:size(subclones)[1]
						for k in 1:size(subclones)[2]
							if (subclones[j,k] == subclones[coordX+x , coordY+y])
								subclones[j,k] = subclones[coordX , coordY]

							end
						end
					end

				end
			end
		end
	end

	if !(assigned)
		global new += 1
		subclones[coordX , coordY] = new
	end

end


# Clean up
for q in 1:size(subclones)[1]
	for r in 1:size(subclones)[2]
		
		for x in -1:1
			for y in -1:1

				if ((x == 0) && (y == 0)) continue end

				# Take care of boundaries
				if ( (q+x < 1) || (q+x > sizeX) || (r+y < 1) || (r+y > sizeY) )
					continue
				end

				if ( (subclones[q,r] != subclones[q+x , r+y]) && (subclones[q+x , r+y] != 0) && (subclones[q,r] != 0) )
					for j in 1:size(subclones)[1]
						for k in 1:size(subclones)[2]
							if (subclones[j,k] == subclones[q+x , r+y])
								subclones[j,k] = subclones[q , r]

							end
						end
					end
				end

			end
		end

	end
end


# Tidy up sub clone labels
#sort_labels(sublcones)



print("\rSeparating sub-clones... done.")



#================= Find boundaries which lie within larger closed boundaries ===================#


println("")
print("Finding internal boundaries...")


for inner_label in 1:findmax(subclones)[1]


	cell_i = 0
	cell_j = 0


	# Find cell in sublcones with this label
	for i in 1:size(subclones)[1]
		for j in 1:size(subclones)[2]

			if (subclones[i,j] == inner_label)
				cell_i = i
				cell_j = j
			end
		end
	end


	# If labels have been merged/deleted cell_i and cell_j will be zero. Continue to next label...
	if ((cell_i == 0) && (cell_j == 0)) continue end
	

	# For this subclone boundary, ask if it is bound to by cells 
	# from any other boundary in the upwards, downwards, left and right direction
	candidates_u = []
	candidates_d = []
	candidates_l = []
	candidates_r = []
	for i in 1:size(subclones)[1]

		if (i < cell_i)
			if ((subclones[i,cell_j] != 0) && (subclones[i,cell_j] != inner_label))
				push!(candidates_l , subclones[i,cell_j])
			end
		end

		if (i > cell_i)
			if ((subclones[i,cell_j] != 0) && (subclones[i,cell_j] != inner_label))
				push!(candidates_r , subclones[i,cell_j])
			end
		end
	end

	for j in 1:size(subclones)[2]

		if (j < cell_j)
			if ((subclones[cell_i,j] != 0) && (subclones[cell_i,j] != inner_label))
				push!(candidates_d , subclones[cell_i,j])
			end
		end

		if (j > cell_j)
			if ((subclones[cell_i,j] != 0) && (subclones[cell_i,j] != inner_label))
				push!(candidates_u , subclones[cell_i,j])
			end
		end
	end


	# Then check if a label appears in all of the candidates_u, candidates_d, candidates_l, candidates_r lists
	outer_label = 0
	for label_j in 1:findmax(subclones)[1]

		if ((label_j in candidates_u) && (label_j in candidates_d) && (label_j in candidates_l) && (label_j in candidates_r))
			outer_label = label_j
		end
	end



	# If no labels are in all 4 lists, move onto next subclone
	if (outer_label == 0) continue end



	# If a label is common to all 4 lists, record path from outer surface of our subclone to this potential "outer" subclone
	chord = []
	for i in cell_i:size(subclones)[1]

			if (subclones[i,cell_j] == inner_label) chord = []
			else push!(chord , sliced_tumour[i+minX-1,cell_j+minY-1]) end

			if (subclones[i,cell_j] == outer_label) break end

	end



	# Then check that the region between inner boundary and outer boundary is red
	if (chord[1] == 1.0)


		# If so, either delete inner boundary or assign it the same label as outer boundary
		if (death == '1')

			# delete internal boundary
			for i in 1:size(subclones)[1]
				for j in 1:size(subclones)[2]

					if (subclones[i,j] == inner_label) subclones[i,j] = 0.0 end

				end
			end

		else

			# assign the same label to inner and outer boundary
			merged_label = min(inner_label , outer_label)
			binned_label = max(inner_label , outer_label)

			for i in 1:size(subclones)[1]
				for j in 1:size(subclones)[2]

					if (subclones[i,j] == binned_label) subclones[i,j] = merged_label end

				end
			end

		end
	end

end


# Tidy up sub clone labels
sort_labels(subclones)


print("\rFinding internal boundaries... Done.")


#================== Write sub-clone boundary data ===================#
println("")
print("Writing data...")


# Write 2D sliced data
outfile = string(path , "/" , slice , "_raw.csv" )
open(outfile, "w") do f
	for i in 1:size(sliced_tumour)[1]
		for j in 1:size(sliced_tumour)[2]
			if ( sliced_tumour[i,j] != -1.0 ) write( f , "$i,$j,$(sliced_tumour[i,j])\n") end
		end
	end
end


# Write boundary data
outfile = string(path , "/" , ARGS[2] , "_sepBoundaries.csv" )
open(outfile, "w") do f

	for i in 1:size(subclones)[1]
		for j in 1:size(subclones)[2]

			if (subclones[i,j] > 0.0)
				write( f , "$(i+minX-1),$(j+minY-1),$(subclones[i,j])\n")
			end
		end
	end
end

print("\rWriting data... done.")
println("")
println("")











