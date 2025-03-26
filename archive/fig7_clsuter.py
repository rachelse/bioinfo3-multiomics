import sys
import numpy as np

# Read in tsv file
# First line is the header
# Fill in matrix (2, 105)
mat = np.zeros((2, 105))

with open(sys.argv[1]) as f:
    line = f.readline()
    names = line.strip().split("\t")
    line = f.readline()
    line = line.strip().split("\t")
    for i in range(1, len(line)):
        mat[0, i-1] = float(line[i])
    line = f.readline()
    line = line.strip().split("\t")
    for i in range(1, len(line)):
        mat[1, i-1] = float(line[i])

# Get argmax of each column
member = np.argmax(mat, axis=0) + 1
membership = mat / np.sum(mat, axis=0)

# Write out the result
with open(sys.argv[2], "w") as f:
    f.write("Cluster\tName\tmembership1\tmembership2\n")
    for i in range(len(member)):
        f.write(f"{names[i]}\t{member[i]}\t{membership[0,i]}\t{membership[1,i]}\n")