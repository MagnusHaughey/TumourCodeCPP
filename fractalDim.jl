


# Round and parse a given floating point to integer
@everywhere function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end




# Compute the minimum number of boxes of size=grid needed to cover the pixels at cooredinates contained within Xinput and Yinput
@everywhere function compute_boxCount(Xinput , Yinput , grid)


	min_boxCount = (length(Xinput))^2


	# Repeat box counting procedure for many different offsets in both x and y
	for xOffset in 0:grid

		for yOffset in 0:grid

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


			# Sum the number of bins with a non-zero pixel count
			boxCount = 0
			for i in 1:size(boxweights)[1]
				for j in 1:size(boxweights)[2]
					if (boxweights[i,j] > 0) boxCount += 1
					elseif (boxweights[i,j] > grid^2) println("Problem with boxcount - check bins. Exiting..."); exit(0) end
			    end
			end

			# Check if we have a new minimum box count for this x-y offset
			if (boxCount < min_boxCount) min_boxCount = boxCount end

		end

	end


	return min_boxCount


end




@everywhere function fractalData(Xinput , Yinput , grids , outpath , label )

	println("")
	boxCount_futures = Array{Future}(undef , length(grids))
	boxCount = zeros(length(grids))

	index = 1
	core = 1
	for grid in grids

		# Compute boxcounts for each gridsize
		boxCount_futures[index] = remotecall( compute_boxCount , core , Xinput , Yinput , grid )

		index += 1
		core %= nworkers()
		core += 1

	end

	# Fetch the values of all boxcounts
	for i in 1:length(grids)

		outval = fetch(wait(boxCount_futures[i]))
		println("Grid=$(grids[i]); boxcount=$outval")
		boxCount[i] = outval

	end


	outfile = string(outpath , "/subclone#$label" , ".dat" )
	open(outfile, "w") do f
		for i in 1:length(grids)
			write( f , "$(-log(grids[i])) $(log(boxCount[i]))\n" )
		end
	end


	# Calcualte fractal dimension
	grad , yint , r , p , stdErr = st.linregress( np.log(grids) , np.log(boxCount) )


	return -grad

end



#= Calculate fractal dimension data point for a given grid size m
*** Note that the (n+1)th index in the first axis of the Qpr array 
corresponds to a boxweight, p, of n. This is to avoid indexing the 
zeroth element of Qpr for boxweight = 0
=#


@everywhere function compute_marginalQpr( input , subclone_size , xBound , yBound , grid )

	marg_Qpr = fill(0.0 , (subclone_size+1))

	Nr = 0
	xOffset = 0

	while( (grid + xOffset - 1) < xBound )

		yOffset = 0

		while( (grid + yOffset - 1) < yBound )
			boxweight = 0
				
			xBegin = xOffset
			yBegin = yOffset
			xEnd = grid + xOffset - 1
			yEnd = grid + yOffset - 1

			for i in xBegin:xEnd
				for j in yBegin:yEnd
					if ([i,j] in input) boxweight += 1 end
				end
			end

			Nr += 1
			marg_Qpr[boxweight+1] += 1

			yOffset += 1
		end

		xOffset += 1
	end

	
	#print("GRID=$grid\n")
	# Normalise distribution
	for p in 1:size(marg_Qpr)[1]
		marg_Qpr[p] /= Nr
	end

	#print(marg_Qpr)

	return marg_Qpr



end


@everywhere function lacunarityData( input , subclone_size , xBound , yBound , grids )


	Qpr_futures = Array{Future}(undef , length(grids))
	#for i in 1:length(Qpr_futures)
	#	Qpr_futures[i] = Array{Future}(undef , (subclone_size+1))
	#end

	Qpr = fill(0.0 , (subclone_size+1 , length(grids)))


	index = 1
	core = 1
	for grid in grids

		# Compute boxcounts for each gridsize (*** Notice that the Qpr_futures array is the "transpose" of the Qpr array)
		Qpr_futures[index] = remotecall( compute_marginalQpr , core , input , subclone_size , xBound , yBound , grid )

		index += 1
		core %= nworkers()
		core += 1

	end

	# Fetch all values of Qpr array when jobs have finished
	for r in 1:length(grids)

		marg = fetch( wait(Qpr_futures[r]) )
		#print(marg)
		#exit(0)

		for p in 1:(subclone_size+1)

		val = marg[p,1]
		Qpr[p , r] = val

		end
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
print("Importing Python/NumPy modules [$(nworkers()) thread(s)]...")
@everywhere using PyCall
@everywhere @pyimport numpy as np
@everywhere @pyimport scipy.stats as st

print("\rImporting Python/NumPy modules [$(nworkers()) thread(s)]... done.")




#================== Read in import data file from command line ===================#
println("")
print("Reading tumour data...")
path = ARGS[1]

# Create directory for fractal dimension and lacunarity data files
frac_path = string( path , "/" , ARGS[2] , "_fractalDimensions/")
lac_path = string( path , "/" , ARGS[2] , "_lacunarity/")

# If directories do not exist create one
if (!isdir(frac_path)) mkdir(frac_path) end
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


	# Construct new array for individual sub-clone
	for i in 1:size(subclones)[1]
		for j in 1:size(subclones)[2]

			if (subclones[i,j] == label) 
				push!(sub , [i,j]) 
				push!(subcloneX , i)
				push!(subcloneY , j)

				if (i < minX) minX = i
				elseif (i > maxX) maxX = i
				end

				if (j < minY) minY = j
				elseif (j > maxY) maxY = j
				end
			end

		end
	end

	if (length(sub) < 9) continue end



	# Re-normalise coordinates of isolated sub-clone
	for gland in sub
		gland[1] -= minX
		gland[2] -= minY
	end

	for i in 1:length(subcloneX)

		subcloneX[i] -= (minX - 1)
		subcloneY[i] -= (minY - 1)

	end

	maxX -= minX
	maxY -= minY
	minX = 0
	minY = 0



	# Compute fractal dimension of given sub-clone
	maxGrid = min(maxX , maxY)

	

	#======= FOR SIERPINSKI CARPET ONLY (AT THE MOMEMT) ========

	# Set list of grid sizes as the scales corresponding to the sequence of minima...
	iteration = 0
	while((3^iteration) < (length(subcloneX)^0.5))
		push!(grids, 3^iteration)
		iteration += 1
	end
	#println(grids)
	#exit(0)

	===========================================================#



	#======= FOR ALL OTHER STRUCTURES ========#

	iteration = 0
	while((2^iteration) < maxGrid)
		push!(grids, 2^iteration)
		iteration += 1
	end

	#=========================================#

#=
	# Construct list of box sizes to be used in fractal dimension calculation
	for i in 1:maxGrid
		push!(grids, i)
	end
=#

	# Compute fractal dimension
	fdim = fractalData(subcloneX , subcloneY , grids , frac_path , label)

	# Append fractal dimension value to array
	push!( fractalDims , ( label , fdim ) )


	println("")
	print("\rCalculating lacunarity curves...")


	# Compute frequency distribution for the number of boxes of size r with a box-mass of p
	Qpr = lacunarityData( sub, length(sub) , maxX , maxY , grids )


	# Compute moments of Qpr distribution and write lacunarity data to file
	outfile = string(lac_path , "/subclone#$label" , ".dat" )
	open(outfile, "w") do f
		for r in 1:length(grids)
			#if (r%20 != 0) continue end
			# Compute first moment of Qpr
			Z1 = pMoments( Qpr , r , 1 )
			# Compute second moment of Qpr
			Z2 = pMoments( Qpr , r , 2 )
			lac = Z2/(Z1^2)
			write( f , "$(log(grids[r])) $(log(lac))\n" )
		end
	end


	print("\rCalculating lacunarity curves... Done.")


end

# Write fractal dimension data
outfile = string(frac_path , "/fractalDimensions.csv" )
open(outfile, "w") do f

	for subclone in fractalDims
		write( f , "$(subclone[1]),$(subclone[2])\n")
	end
end



println("\n")



















