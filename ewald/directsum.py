import numpy as np

r = np.zeros([3, 8])
L = 10.0
q = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]

#r[0][4] = 0.25 * L
#r[1][4] = 0.25 * L
#r[2][4] = 0.25 * L
#r[0][5] = 0.25 * L
#r[1][5] = 0.75 * L
#r[2][5] = 0.75 * L
#r[0][6] = 0.75 * L
#r[1][6] = 0.25 * L
#r[2][6] = 0.75 * L
#r[0][7] = 0.75 * L
#r[1][7] = 0.75 * L
#r[2][7] = 0.25 * L

r[0][4] = 0.5 * L
r[1][5] = 0.5 * L
r[2][6] = 0.5 * L
r[0][7] = r[1][7] = r[2][7] = 0.5 * L

r[0][1] = 0.5 * L
r[1][1] = 0.5 * L
r[0][2] = 0.5 * L
r[2][2] = 0.5 * L
r[1][3] = 0.5 * L
r[2][3] = 0.5 * L

MAX_N = int(raw_input('Enter MAX_N: '))

energy = 0.0

rij = np.zeros(3)
for i in range(8):
    for j in range(8):
        if ( i != j ):
            for k in range(3):
                rij[k] = r[k][i] - r[k][j]
            drij = np.sqrt(np.dot(rij, rij))
            energy += q[i]*q[j]/drij
                
rn = np.zeros(3)
for ni in range(-MAX_N, MAX_N):
    for nj in range(-MAX_N, MAX_N):
#        for nk in range(-MAX_N, MAX_N):
            
        rn[0] = ni*L
        rn[1] = nj*L
        rn[2] = 0.0
        if (np.dot(rn, rn) != 0):
            
            for i in range(8):
                for j in range(8):
                    for k in range(3):
                        rij[k] = r[k][j] - r[k][i] + rn[k]
                    drij = np.sqrt( np.dot(rij, rij) )
                    energy += q[i]*q[j]/drij

energy = 0.5 * energy
                    
                                
                

print energy
            
            
