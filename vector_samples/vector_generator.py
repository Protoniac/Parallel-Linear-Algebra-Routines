import numpy as np
import sys

size=int(sys.argv[1])

f = open(str(size)+"_"+sys.argv[2]+".txt","w")
vector = np.random.rand(size)*10
for i in range(size):
    f.write((str(round(vector[i],4))+"\n"))
