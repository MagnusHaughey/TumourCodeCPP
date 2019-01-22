



# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end


#=
# Calculate fractal dimension data point for a given grid size m
function fractalData( input , xBound , yBound , m )

	xBoxes = 0
	yBoxes = 0

	N_orientations = 0
	N_cumulative = 0
	N = 0

	yOffset = 0

		while(yOffset <= m)

			print("\rBox size = $m; y-offset = $yOffset")

			N_orientations += 1
			N = 0

			while( (xBoxes*m) < xBound )

				while( (yBoxes*m+yOffset) <= yBound )


					if (yBoxes == 0)

						gland_found = false

						for i in (xBoxes*m):(m-1+(xBoxes*m))
							if gland_found break end

							for j in 0:(yOffset-1)

								if gland_found break end

								if [i,j] in input
									N += 1
									#println("Gland found!\n")
									gland_found = true
								end
							end
						end
					end

					

						gland_found = false

						for i in (xBoxes*m):(m-1+(xBoxes*m))
							if gland_found break end

							for j in (yBoxes*m+yOffset):(m-1+(yBoxes*m)+yOffset)

								if gland_found break end

								if [i,j] in input
									N += 1
									#println("Gland found!\n")
									gland_found = true
								end
							end
						end
					


				yBoxes += 1

				end

				xBoxes += 1
				yBoxes = 0

			end

			N_cumulative += N
			xBoxes = 0
			yBoxes = 0
			yOffset += 1

		end


	return (N_cumulative/N_orientations)

end
=#


function fractalData2(Xinput , Yinput , grids )

	println("")
	Nboxes = []

	for grid in grids

		print("\rBox size = $grid of $(np.max(grids)) ")

		N_orientations = 0
		N_cumulative = 0

		min_boxcount = (length(Xinput))^2

		# Repeat box counting procedure for many different offsets in both x and y
		for xOffset in 0:grid

			for yOffset in 0:grid

				N_orientations += 1

				Xinput_offsetted = []
				Yinput_offsetted = []

				for k in 1:length(Xinput)
					push!(Xinput_offsetted , (Xinput[k] + xOffset))
					push!(Yinput_offsetted , (Yinput[k] + yOffset))
				end

				# Compute 2 dimensional histogram of subclone boundary data 
				xBins = np.arange(0,(np.max(Xinput)[1]+xOffset+3+grid),grid)
				yBins = np.arange(0,(np.max(Yinput)[1]+yOffset+3+grid),grid)
				boxweights, xEdges , yEdges = np.histogram2d( Xinput_offsetted , Yinput_offsetted , bins=(xBins,yBins) )


#=
				if (grid >= 60)

					#println("$(size(boxweights))")
					tot = 0
					for i in 1:size(boxweights)[1]
						for j in 1:size(boxweights)[2]
							println("[$i,$j] -> $(boxweights[i,j])")
							tot += boxweights[i,j]
					    end
					end
					if (tot != length(Xinput)) println("Missing some pixels! Exiting...\n"); exit(0) end
					println("\n\n")
					#print(xEdges)
					#println("\n\n")
					#print(yEdges)
					#println("\n\n")
					#println("$(np.max(Xinput)[1]) $(np.max(Yinput)[1])")
					#exit(0)
				end
=#


				# Sum the number of bins with a non-zero pixel count
				boxcount = 0
				for i in 1:size(boxweights)[1]
					for j in 1:size(boxweights)[2]
						if (boxweights[i,j] > 0) boxcount += 1
						elseif (boxweights[i,j] > grid^2) println("Problem with boxcount - check bins. Exiting..."); exit(0) end
				    end
				end

				# Check if we have a new minimum box count for this x-y offset
				if (boxcount < min_boxcount) min_boxcount = boxcount end

			end

		end

		#println("Grid = $grid -> minimum count = $min_boxcount\n")

		push!(Nboxes , min_boxcount)
	    
	end


	# Calcualte fractal dimension
	grad , yint , r , p , stdErr = st.linregress( np.log(grids) , np.log(Nboxes) )


	return -grad


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





#================== Compute fractal dimension of each sub-clone ===================#
println("")
print("Calculating fractal dimensions...")

fractalDims = []

for label in 1:findmax(subclones)[1]
	sub = []
	dfData = []
	lacData = []
	grids = []
	subcloneX = []
	subcloneY = []
	minX = size(subclones)[1]
	minY = size(subclones)[2]
	maxX = 0
	maxY = 0


	# Construct new array for individual sub-clone (two lists with x- and y-coordinates)
	for i in 1:size(subclones)[1]
		for j in 1:size(subclones)[2]

			if (subclones[i,j] == label) 
				push!(subcloneX , i)
				push!(subcloneY , j)

				# Update max and min x- and y-coordinates
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


	# Re-normalise coordinates of isolated sub-clone
	for gland in sub
		gland[1] -= minX
		gland[2] -= minY
	end
	maxX -= minX
	maxY -= minY
	minX = 0
	minY = 0


	# Define maximum gridsize for fractal dimension calcualtions
	maxGrid = min(maxX , maxY)


	# Construct list of box sizes to be used in fractal dimension calculation
	for i in 1:maxGrid
		push!(grids, i)
	end


	# Compute fractal dimension
	fdim = fractalData2(subcloneX , subcloneY , grids)

	println("")


	# Compute frequency distribution for the number of boxes of size r with a box-mass of p
	#Qpr = lacunarityData( sub, length(sub) , maxX , maxY , maxGrid )


	# Append fractal dimension value to array
	push!( fractalDims , ( label , fdim ) )

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



















