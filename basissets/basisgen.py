BASISNAME = "standard.basis"
BASE = 2.0
MINEXP = -3
NEXPS = 7

outfile = open(BASISNAME, 'w')

nbfs = 0
alpha = 0.0
for a in range (NEXPS+1):
    beta = 0.0
    for b in range (NEXPS+1):
        gamma = 0.0
        for c in range (NEXPS+1):
            outfile.write("basisfunction\nalpha,\n")
            outfile.write("0.0, %f\n0.0, 0.0\nbeta,\n" % (alpha) )
            outfile.write("0.0, 0.0\n%f, 0.0\ngamma,\n" %(beta) )
            outfile.write("0.0, %f\n0.0, 0.0\nendfunction\n\n" %(gamma))
            nbfs = nbfs + 1
            gamma = BASE**(MINEXP + c)
        beta = BASE**(MINEXP + b)
    alpha = BASE**(MINEXP + a)

print nbfs


