import numpy as np
import sys

nrow=int(sys.argv[1])
ncol=int(sys.argv[2])

f = open(str(nrow)+"x"+str(ncol)+"_"+sys.argv[3]+".txt","w")
matrix = np.random.rand(nrow,ncol)*10
for i in range(nrow):
    for j in range(ncol):
        if(j!=ncol-1):
            f.write((str(round(matrix[i,j],4))+",".rstrip('\n')))
        else:
            f.write((str(round(matrix[i,j],4))+"".rstrip('\n')))
    f.write('\n')
