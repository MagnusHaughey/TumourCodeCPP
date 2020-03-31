

import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import sys


infile = sys.argv[1]

# Import raw data
x , y , mutations = np.loadtxt(infile , delimiter="," , unpack=True)
#x_boundary , y_boundary , mutations_boundary = np.loadtxt(infile + '.sepBoundaries.csv' , delimiter="," , unpack=True)




maxX = int(np.max(x))
maxY = int(np.max(y))

minX = int(np.min(x))
minY = int(np.min(y))



H = np.ones( (maxX-minX , maxY-minY) )
H = -H

# Fill array H
for i in range(len(x)):
	if ( int(mutations[i]) >= 1 ):
		H[int(x[i])-minX-1 , int(y[i])-minY-1] = 1

	else: 
		H[int(x[i])-minX-1 , int(y[i])-minY-1] = int(mutations[i])




#for i in range(len(x_boundary)):
#	H[int(x_boundary[i])-minX-1 , int(y_boundary[i])-minY-1] = 2




# Adjust bins so that final image is always square
rangeX = maxX - minX
rangeY = maxY - minY

#print(rangeX)
#print(rangeY)

if (rangeX >= rangeY):
	xStart = minX
	xEnd = maxX
	xPartitions = rangeX + 1

	yStart = minY + 0.5*(rangeY - rangeX)
	yEnd = minY + 0.5*(rangeY + rangeX)
	yPartitions = rangeX + 1

else:
	yStart = minY
	yEnd = maxY
	yPartitions = rangeY + 1

	xStart = minX + 0.5*(rangeX - rangeY)
	xEnd = minX + 0.5*(rangeX + rangeY)
	xPartitions = rangeY + 1

'''
xStart = minX
xEnd = maxX
xPartitions = maxX - minX + 1
yStart = minY
yEnd = maxY
yPartitions = maxY - minY + 1
'''

xedges = np.linspace(xStart , xEnd , xPartitions)
yedges = np.linspace(yStart , yEnd , yPartitions)


H = H.T



if (len( [val for val in mutations if (val == 1)] ) > 0):
	cmap = mpl.colors.ListedColormap(['white', 'green', 'red'])
	#cmap = mpl.colors.ListedColormap(['white', 'green', 'red', 'black'])
else:
	cmap = mpl.colors.ListedColormap(['white', 'green'])

#plt.switch_backend('agg')

plt.imshow(H, interpolation='none' , origin='lower' , cmap=cmap , extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.axis('off')
#plt.show()
plt.savefig(infile + ".png" , dpi=1000 , format='png')
#print(H)











