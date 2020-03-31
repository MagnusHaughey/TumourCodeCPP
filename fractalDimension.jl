


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
					elseif (boxweights[i,j] > grid^2) #println("Problem with boxcount - check bins. Exiting...");
						println(0)
						exit(0)
					end
			    end
			end

			# Check if we have a new minimum box count for this x-y offset
			if (boxCount < min_boxCount) min_boxCount = boxCount end

		end

	end


	return min_boxCount


end




@everywhere function fractalData(Xinput , Yinput , grids , outpath , label )

	#println("")
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
		#println("Grid=$(grids[i]); boxcount=$outval")
		boxCount[i] = outval

	end


	outfile = string(outpath , "/subclone#$label" , ".dat" )
	open(outfile, "w") do f
		for i in 1:length(grids)
			write( f , "$(-log(grids[i])) $(log(boxCount[i]))\n" )
		end
	end

	#=
	# Calcualte fractal dimension
	outfile = string(frac_path , "../../errors.dat" )
	open(outfile, "a") do f

		writedlm(f , "$seed: problem with linear regression. length(grids)=$(length(grids)), length(boxCount)=$(length(boxCount))\n")

	end
	=#
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




#=====================================================================#



# Determine quiet or verbose
if (length(ARGS) < 2)
	quiet = true
elseif string(ARGS[1]) == "-q"
	quiet = true
else
	quiet = false
end


#================== Import Python/NumPy modules ===================#
if (!quiet) println("") end 
if (!quiet) print("Importing Python/NumPy modules [$(nworkers()) thread(s)]...") end
@everywhere using PyCall
#@everywhere using CSV
@everywhere using DelimitedFiles
@everywhere @pyimport numpy as np
@everywhere @pyimport scipy.stats as st

if (!quiet) print("\rImporting Python/NumPy modules [$(nworkers()) thread(s)]... done.") end




#================== Read in import data file from command line ===================#
if (!quiet) println("") end
if (!quiet) print("Reading tumour data...") end

if (length(ARGS) < 2)
	path = ARGS[1]
else 
	path = ARGS[2]
end


# Interpret simulation parameters
split1 = split(path , "/")
for dirs in split1
	if occursin("maxsize" , dirs)

		split2 = split(dirs , "_")

		for fragment in split2

			if occursin("seed" , fragment)
				global seed = parse(Int64 , split(fragment, "=")[2])
		    end
		end

	end
end


# Create directory for fractal dimension and lacunarity data files
#frac_path = string( path , "/" , slice , "_fractalDimensions/")
frac_path = string( path , ".fractalDimensions/" )
#lac_path = string( path , "/" , ARGS[2] , "_lacunarity/")

# If directories do not exist create one
if (!isdir(frac_path)) mkdir(frac_path) end
#if (!isdir(lac_path)) mkdir(lac_path) end

infile = string(path)
data = np.loadtxt( infile , unpack=true , delimiter="," , usecols=(0,1,2))

# Construct and populate subclones array
minX = findmax(data)[1]
minY = findmax(data)[1]
maxX = 0
maxY = 0

# First line indicates if resistant/WT boundary is the first (1,3,4) or second (2,5,6) shape
#RES_FIRST = false
#if (data[1] == 1)
#	RES_FIRST = true
#end


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
subclones_combined = fill(0.0 , (sizeX,sizeY))





# Unpack into separate x, y, z, GA vectors
for i in 1:length(data)
	if ( ((i-3)%3) == 0 )
		subclones[float_to_int(data[i-2]-minX+1), float_to_int(data[i-1]-minY+1)] = data[i]
	end
end

if (!quiet) print("\rReading tumour data... done.") end





#================== Compute fractal dimension & lacunarity of each sub-clone edge ===================#
if (!quiet) println("") end
if (!quiet) print("Calculating fractal dimensions...") end

fractalDims = []
imageSize_list = []
boundaries_to_process = []




# Compute fractal dimension for each individual boundary
for label in 1:findmax(subclones)[1]

	# Construct new array for individual sub-clone	
	subcloneX = [index[1] for index in findall(subclones .== label)]
	subcloneY = [index[2] for index in findall(subclones .== label)]

	push!(boundaries_to_process , [subcloneX , subcloneY , label])

end







for boundary in boundaries_to_process

	grids = []

	subcloneX = boundary[1]
	subcloneY = boundary[2]
	label = boundary[3]



	if (length(subcloneX) < 8) continue end

	minX = findmin(subcloneX)[1]
	maxX = findmax(subcloneX)[1]

	minY = findmin(subcloneY)[1]
	maxY = findmax(subcloneY)[1]



	for i in 1:length(subcloneX)

		subcloneX[i] -= (minX - 1)
		subcloneY[i] -= (minY - 1)

	end




	#======= Create list of grids ========#

	#maxGrid = float_to_int(0.005*length(subcloneX))
	maxGrid = max(maxX , maxY)
	push!( imageSize_list , maxGrid )

	for i in 2:maxGrid

		push!(grids , i)

	end

	if (length(grids) < 2)			# Make sure there are at least 2 grids used 
		push!(grids , grids[1]+1)
	end



	if (length(subcloneX) == 0)

		println(0)
		exit(0)

	end

	

	# Compute fractal dimension
	fdim = fractalData(subcloneX , subcloneY , grids , frac_path , label)

	# Append fractal dimension value to array
	push!( fractalDims , ( label , fdim ) )




	#println("\nsubclone #$(label): imageSize=$imageSize; arm length = $(length(subcloneX)); grids=$grids; --> Df = $fdim")

#=

	if (!quiet) println("") end
	if (!quiet) print("\rCalculating lacunarity curves...") end




	#======= Create list of grids ========#

	maxGrid = float_to_int(0.5*length(subcloneX))
	push!( imageSize_list , maxGrid )

	#increment = max( 1 , floor(maxGrid/(imageSize_list[1])) )
	increment = 2

	for i in 2:increment:maxGrid

		push!(grids , i)
	end


	# Compute frequency distribution for the number of boxes of size r with a box-mass of p
	Qpr = lacunarityData( sub, length(sub) , maxX , maxY , grids )


	# Compute moments of Qpr distribution and write lacunarity data to file
	outfile = string(lac_path , "/subclone#$label" , ".dat" )
	open(outfile, "w") do f
		for r in 1:length(grids)

			# Compute first moment of Qpr
			Z1 = pMoments( Qpr , r , 1 )
			# Compute second moment of Qpr
			Z2 = pMoments( Qpr , r , 2 )

			lac = Z2/(Z1^2)
			write( f , "$(log(grids[r])) $(log(lac))\n" )
		end
	end


	if (!quiet) print("\rCalculating lacunarity curves... Done.") end
=#

end


if (!quiet) print("\rCalculating fractal dimensions... Done.\n") end




#========================================================================#


# Write fractal dimension data to separate file in case appending data to large file throws an error
outfile = string(frac_path , "/fractalDimensions.csv" )
open(outfile, "w") do f

	for i in 1:length(fractalDims)
		write( f , "$(fractalDims[i][2])\n")
	end

end


#=
# Also write fractal dimension data to a big file with all values for all tumours
outfile = string(frac_path , "../../all_dimensions.dat" )
open(outfile, "a") do f

	to_write = []
	for i in 1:length(fractalDims)
		push!(to_write , fractalDims[i][2])
	end

	writedlm(f , [to_write] , " ")
	#writedlm(f , [[fractalDims[1][2] , fractalDims[2][2] , fractalDims[3][2]]] , " ")

end
=#




#println(1)
















