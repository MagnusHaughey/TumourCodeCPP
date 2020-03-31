

import numpy as np
import pandas as pd
import sys
import math


rawData = pd.read_csv(sys.argv[1] , delimiter=" " , header=None)
rawData = rawData.values


for i in range(int(np.max([line[0] for line in rawData])) + 1):

	all_values = [line[1] for line in rawData if (line[0] == i) and (np.isnan(line[1]) == False)]

	print("{} {} {}".format( i , np.mean(all_values) , np.std(all_values)/math.sqrt(len(all_values)) ))	
