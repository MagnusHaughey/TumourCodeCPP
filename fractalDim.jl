



# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end



# Calculate fractal dimension data point for a given grid size m
function fractalData( input , xBound , yBound , m )

	N = 0

	xBoxes = 0
	yBoxes = 0

	while( (xBoxes*m) < xBound )

		while( (yBoxes*m) < yBound )

		gland_found = false

			for i in (xBoxes*m):(m-1+(xBoxes*m))
				if gland_found break end

				for j in (yBoxes*m):(m-1+(yBoxes*m))

					if gland_found break end

					if [i,j] in input
						N += 1
						gland_found = true
					end

				end

			end

			yBoxes += 1

		end

		xBoxes += 1
		yBoxes = 0

	end

	#println("m=$m , N=$N\n")

	return N

end


#= Calculate fractal dimension data point for a given grid size m

*** Note that the (n+1)th index in the first axis of the Qpr array 
corresponds to a boxweight, p, of n. This is to avoid indexing the 
zeroth element of Qpr for boxweight = 0

=#

function lacunarityData( input , subclone_size , xBound , yBound , maxGrid )


	Qpr = fill(0.0 , (subclone_size+1 , maxGrid))

	for gridsize in 1:maxGrid

		#if (gridsize%20 != 0) continue end

		print("\rCalculating lacunarity... [$gridsize/$maxGrid]")

		Nr = 0
		xOffset = 0

		while( (gridsize + xOffset - 1) < xBound )

			yOffset = 0

			while( (gridsize + yOffset - 1) < yBound )
				boxweight = 0
				
				xBegin = xOffset
				yBegin = yOffset
				xEnd = gridsize + xOffset - 1
				yEnd = gridsize + yOffset - 1

				for i in xBegin:xEnd
					for j in yBegin:yEnd
						if ([i,j] in input) boxweight += 1 end
					end
				end

				Nr += 1
				Qpr[boxweight+1 , gridsize] += 1

				yOffset += 1
			end

			xOffset += 1
		end

		# Normalise distribution
		for p in 1:size(Qpr)[1]
			Qpr[p , gridsize] /= Nr
		end

		#println("r=$gridsize -> N(r)=$Nr")

	end

	

	return Qpr

end



function pMoments( input , fixed_r , moment )


	outval = 0.0

	# Use trapezoid method to compute moment
	for p in 1:size(input)[1]

		outval += ((p-1)^moment) * input[p , fixed_r]

	end


	return outval

end





#================== Import Python/NumPy modules ===================#
println("")
print("Importing Python/NumPy modules...")
using PyCall
@pyimport numpy as np
@pyimport scipy.stats as st

print("\rImporting Python/NumPy modules... done.")




#================== Read in import data file from command line ===================#
println("")
print("Reading tumour data...")
path = ARGS[1]

# Create directory for lacunarity data files
lac_path = string( path , "/" , ARGS[2] , "_lacunarity/")

# If directory does not exist create one
if (!isdir(lac_path)) mkdir(lac_path) end

infile = string(path , "/" , ARGS[2] , "_sepBoundaries.csv")
data = np.loadtxt( infile , unpack=true , delimiter="," , usecols=(0,1,2))

# Construct and populate subclones array
minX = findmax(data)[1]
minY = findmax(data)[1]
maxX = 0
maxY = 0

for i in 1:length(data)
	if ( ((i-3)%3) == 0 )
		if (data[i-2] < minX) global minX = data[i-2] end
		if (data[i-2] > maxX) global maxX = data[i-2] end
		if (data[i-1] < minY) global minY = data[i-1] end
		if (data[i-1] > maxY) global maxY = data[i-1] end
	end
end


sizeX = float_to_int(maxX - minX) + 1
sizeY = float_to_int(maxY - minY) + 1
subclones = fill(0.0 , (sizeX,sizeY))

# Unpack into separate x, y, z, GA vectors
for i in 1:length(data)
	if ( ((i-3)%3) == 0 )
		subclones[float_to_int(data[i-2]-minX+1), float_to_int(data[i-1]-minY+1)] = data[i]
	end
end

print("\rReading tumour data... done.")


#=
#================== Separate different sub-clones ===================#
println("")
print("Separating sub-clones...")

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



# Re-sort sub clone labels
for i in 1:(new-2)

	next_label = 0

	if !(i in subclones)
		#println("\n!!! i=$i")

		for q in (i+1):(new-1)
			if (q in subclones) next_label = q end 
		end
		#print("\n$next_label")

		if (next_label != 0)
			for j in 1:size(subclones)[1]
				for k in 1:size(subclones)[2]
					if (subclones[j,k] == next_label) subclones[j,k] = i end
				end
			end
		end
	end
end




print("\rSeparating sub-clones... done.")






#================== Write sub-clone boundary data ===================#
println("")
print("Writing data...")

# Write 2D sliced data
outfile = string(path , "/" , ARGS[2] , "_sep_boundaries.csv" )
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
=#



#================== Compute fractal dimension of each sub-clone ===================#
println("")
print("Calculating fractal dimensions...")

fractalDims = []

for label in 1:findmax(subclones)[1]
	sub = []
	dfData = []
	lacData = []
	grids = []
	minX = size(subclones)[1]
	minY = size(subclones)[2]
	maxX = 0
	maxY = 0

	#println("\n minX=$minX --- minY=$minY --- maxX=$maxX --- maxY=$maxY \n")

	# Construct new array for individual sub-clone
	for i in 1:size(subclones)[1]
		for j in 1:size(subclones)[2]

			if (subclones[i,j] == label) 
				push!(sub , [i,j]) 

				if (i < minX) minX = i
				elseif (i > maxX) maxX = i
				end

				if (j < minY) minY = j
				elseif (j > maxY) maxY = j
				end
			end

		end
	end

	if (length(sub) < 4) continue end

	#println("\n minX=$minX --- minY=$minY --- maxX=$maxX --- maxY=$maxY \n")

	# Re-normalise coordinates of isolated sub-clone
	for gland in sub
		gland[1] -= minX
		gland[2] -= minY
	end
	maxX -= minX
	maxY -= minY
	minX = 0
	minY = 0

	#println("\n minX=$minX --- minY=$minY --- maxX=$maxX --- maxY=$maxY \n")

	# Compute fractal dimension of given sub-clone
	gridsize = min(maxX , maxY)
	maxGrid = gridsize
	
	println("")

	while(gridsize > 0)
		if (gridsize%2 == 0)
			gridsize -= 1
		 	continue
		 end
		print("\rBox size = $gridsize")
		push!(dfData , fractalData( sub , maxX , maxY , gridsize ) )
		push!(grids , gridsize)
		gridsize -= 1
	end


#=
### For particularly large structures, only use a selection of box sizes
	while(gridsize > (gridsize = 100)))
		print("\rBox size = $gridsize")
		push!(dfData , fractalData( sub , maxX , maxY , gridsize ) )
		push!(grids , gridsize)
		gridsize -= 1
	end
	gridsize = 1
	while(gridsize < 100))
		print("\rBox size = $gridsize")
		push!(dfData , fractalData( sub , maxX , maxY , gridsize ) )
		push!(grids , gridsize)
		gridsize += 1
	end
############
=#


	# Compute frequency distribution for the number of boxes of size r with a box-mass of p
	#Qpr = lacunarityData( sub, length(sub) , maxX , maxY , maxGrid )

	# Calcualte fractal dimension
	grad , yint , r , p , stdErr = st.linregress( np.log(grids) , np.log(dfData) )

	# Append fractal dimension value to array
	push!( fractalDims , ( label , -grad ) )

#=
	# Compute moments of Qpr distribution and write lacunarity data to file
	outfile = string(lac_path , "/subclone#$label" , ".dat" )
	open(outfile, "w") do f

		for r in 1:maxGrid

			#if (r%20 != 0) continue end

			# Compute first moment of Qpr
			Z1 = pMoments( Qpr , r , 1 )

			# Compute second moment of Qpr
			Z2 = pMoments( Qpr , r , 2 )

			lac = Z2/(Z1^2)

			write( f , "$(log(r)) $(log(lac))\n" )


		end

	end
=#

end

# Write fractal dimension data
outfile = string(path , "/" , ARGS[2] , "_fractalAnalysis.csv" )
open(outfile, "w") do f

	for subclone in fractalDims
		write( f , "$(subclone[1]),$(subclone[2])\n")
	end
end



print("\rCalculating fractal dimensions... done.")
println("\n")





















