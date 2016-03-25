import numpy as np
import matplotlib.pyplot as plt

kappa = 0.25
param = 0.6
N = 4

def calc_abc(r, R):
    avec = np.zeros([N, N, 3])
    bvec = np.zeros([N, N, 3])
    cvec = np.zeros([N, N, 3])
    avals = np.zeros([N, N])
    bvals = np.zeros([N, N])
    cvals = np.zeros([N, N])

    Rij = np.zeros(3)
    for i in range (N-1):
        Rij[:] = R[i][:]
        for j in range(i+1, N):
            avec[i][j][:] = Rij[:] - r[i][:]
            bvec[j][i][:] = -Rij[:] - r[j][:]
            cvec[i][j][:] = Rij[:] - r[i][:] + r[j][:]

            avals[i][j] = np.sqrt(np.vdot(avec[i][j][:], avec[i][j][:]))
            bvals[j][i] = np.sqrt(np.vdot(bvec[j][i][:], bvec[j][i][:]))
            cvals[i][j] = np.sqrt(np.vdot(cvec[i][j][:], cvec[i][j][:]))
            Rij[:] = Rij[:] + R[j][:]
            
    return [avals, bvals, cvals]

def Jf(z, theta):
    return theta*z/(theta + z)

def calc_wf(r, R):
    avals, bvals, cvals = calc_abc(r, R)
    wf = 0.0
    for i in range(N-1):
        for j in range(i+1, N):
            wf = wf + 0.5*Jf(cvals[i][j], param)
            wf = wf - Jf(avals[i][j], param) - Jf(bvals[j][i], param)
    wf = kappa*wf
    for i in range(N):
        wf = wf - np.vdot(r[i], r[i])
    
    return np.exp(wf)

def calc_grads(r, R, h):
    grads = np.zeros([N, 3])
    rtemp1 = np.zeros([N, 3])
    rtemp2 = np.zeros([N, 3])
    for k in range(N):
        for i in range(3):
            rtemp1[:][:] = r[:][:]
            rtemp2[:][:] = r[:][:]
            rtemp1[k][i] = rtemp1[k][i] + h
            rtemp2[k][i] = rtemp2[k][i] - h

            wf1 = calc_wf(rtemp1, R)
            wf2 = calc_wf(rtemp2, R)

            grads[k][i] = (wf1-wf2)/(2.0*h)
            
    return grads

def calc_laplacian(r, R, h):
    laplacian = 0.0
    rtemp1 = np.zeros([N, 3])
    rtemp2 = np.zeros([N, 3])

    wf = calc_wf(r, R)
    for k in range(N):
        for i in range(3):
            rtemp1[:][:] = r[:][:]
            rtemp2[:][:] = r[:][:]
            rtemp1[k][i] = rtemp1[k][i] + h
            rtemp2[k][i] = rtemp2[k][i] - h
            wf1 = calc_wf(rtemp1, R)
            wf2 = calc_wf(rtemp2, R)
            laplacian = laplacian + wf1 + wf2

    laplacian = (laplacian - 6.0*N*wf)/(h**2)
    return laplacian
        
r = np.zeros([N, 3])
R = np.zeros([N, 3])
r[0][0] = r[0][1] = r[0][2] = 0.5
r[1][0] = -0.1
r[1][1] = 1.0
r[1][2] = 0.5
r[2][0] = 0.4
r[2][1] = 0.5
r[2][2] = -0.5
r[3][0] = 0.5
r[3][1] = -0.1
r[3][2] = 0.6
R[0][0] = -0.2
R[0][1] = 1.0
R[0][2] = 0.8
R[1][0] = 0.7
R[1][1] = -0.3
R[1][2] = -0.6
R[2][0] = 1.0
R[2][1] = 0.0
R[2][2] = 1.0
R[3][0] = R[3][1] = R[3][2] = 0.0


hvals = [0.1, 0.05, 0.001, 0.005, 0.0001, 0.000075, 0.00005, 0.000025, 0.00001, 0.0000075, 0.000005]
lapxs = np.zeros(len(hvals))
gradxs = np.zeros([len(hvals), N, 3])
for i in range(len(hvals)):
    laplacian = calc_laplacian(r, R, hvals[i])
    lapxs[i] = laplacian
    gradxs[i] = calc_grads(r, R, hvals[i])
    
print gradxs
print lapxs

plt.plot(-np.log(hvals), lapxs, 'bd')
final_val = [lapxs[len(hvals)-1]]*len(hvals)
plt.plot(-np.log(hvals), final_val, 'k-')
plt.xlabel('-log h')
plt.ylabel('Laplacian by central differences')
plt.show()




