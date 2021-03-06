# Parameters
PREFIX = "scan2/scan2"

OMEGA = 0.5
MU = 0.5

QMIN = 0.5
QMAX = 2.0
NQ = 31

RMIN = 0.4
RMAX = 10.0
NR = 25

# Main routine
import numpy as np
import matplotlib.pyplot as plt
import sys
qlist = np.linspace(QMIN, QMAX, NQ)
rlist = np.linspace(RMIN, RMAX, NR)

def writeFile(filename, q, r):
    outfile = open(filename, 'w')
    output = "N, 2\nmu, %f, %f\nomega, %f, %f\nq, %f, %f\ngeom,\n%f, 0.0, 0.0\n" % (MU, MU, OMEGA, OMEGA, q, q, r)
    outfile.write(output)
    outfile.close()
    
def extract(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    retval = 0.0
    for line in lines:
        words = line.split()
        if(len(words) == 4):
            retval = words[3]
            break
    infile.close()
    return retval
            
args = sys.argv
s = -1
if (len(sys.argv) == 1):
    print "Please supply an argument"
elif (args[1] == "generate"):
    s = 0
elif (args[1] == "extract"):
    s = 1

for i in range(len(qlist)):
    yvals = np.zeros(len(rlist))
    for j in range(len(rlist)):
        filename = "%s_q%i_r%i" % (PREFIX, i, j)
        if s == 0: # Generate
            filename += ".input"
            writeFile(filename, qlist[i], rlist[j])
        elif s == 1: # Extract
            filename += ".output"
            yvals[j] = extract(filename)
            print qlist[i], rlist[j], yvals[j]
    if s == 1:
        plt.plot(rlist, yvals)

if s == 1:
    plt.show()
            
