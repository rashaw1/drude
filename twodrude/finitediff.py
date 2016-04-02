import numpy as np
import matplotlib.pyplot as plt

muq2 = 2.0*(1.5**2)

def calc_abc(r1, r2, R):
    avec = R - r1
    bvec = -R - r2
    cvec = R - r1 + r2
    a = np.sqrt(np.vdot(avec, avec))
    b = np.sqrt(np.vdot(bvec, bvec))
    c = np.sqrt(np.vdot(cvec, cvec))
    return [a, b, c]

def calc_wf(r1, r2, R):
    a, b, c = calc_abc(r1, r2, R)
    wf = np.exp(-(np.vdot(r1, r1) + np.vdot(r2, r2)))
    wf = wf*np.exp(muq2*0.5*c/(1.0 + c))*np.exp(-muq2*a/(1.0+a))
    wf = wf*np.exp(-muq2*b/(1.0+b))
    return wf

def calc_grad1(r1, r2, R, dt):
    grad1 = np.zeros(3)
    r1temp1 = np.zeros(3)
    r1temp2 = np.zeros(3)
    for i in range(3):
        r1temp1[:] = r1[:]
        r1temp2[:] = r1[:]
        r1temp1[i] = r1temp1[i] + dt
        r1temp2[i] = r1temp2[i] - dt
        wf1 = calc_wf(r1temp1, r2, R)
        wf2 = calc_wf(r1temp2, r2, R)
        grad1[i] = (wf1 - wf2)/(2.0*dt)
    return grad1

def calc_grad2(r1, r2, R, dt):
    grad2 = np.zeros(3)
    r2temp1 = np.zeros(3)
    r2temp2 = np.zeros(3)
    for i in range(3):
        r2temp1[:] = r2[:]
        r2temp2[:] = r2[:]
        r2temp1[i] = r2temp1[i] + dt
        r2temp2[i] = r2temp2[i] - dt
        wf1 = calc_wf(r1, r2temp1, R)
        wf2 = calc_wf(r1, r2temp2, R)
        grad2[i] = (wf1 - wf2)/(2.0*dt)
    return grad2

def calc_laplacian(r1, r2, R, dt):
    laplacian = 0.0
    rtemp1 = np.zeros(3)
    rtemp2 = np.zeros(3)

    wf = calc_wf(r1, r2, R)
    for i in range(3):
        rtemp1[:] = r1[:]
        rtemp2[:] = r1[:]
        rtemp1[i] = rtemp1[i] + dt
        rtemp2[i] = rtemp2[i] - dt
        wf1 = calc_wf(rtemp1, r2, R)
        wf2 = calc_wf(rtemp2, r2, R)
        laplacian = laplacian + wf1 + wf2

        rtemp1[:] = r2[:]
        rtemp2[:] = r2[:]
        rtemp1[i] = rtemp1[i] + dt
        rtemp2[i] = rtemp2[i] - dt
        wf1 = calc_wf(r1, rtemp1, R)
        wf2 = calc_wf(r1, rtemp2, R)
        laplacian = laplacian + wf1 + wf2

    laplacian = (laplacian - 12.0*wf)/(dt**2)
    return laplacian
        
r1 = np.zeros(3)
r2 = np.zeros(3)
R = np.zeros(3)
r1[0] = 0.5
r1[1] = -0.1
r1[2] = 1.0
r2[0] = 0.5
r2[1] = 0.5
r2[2] = -0.5
R[0] = 1.0
R[1] = 1.0
R[2] = 0.0

dtvals = [0.1, 0.05, 0.001, 0.005, 0.0001, 0.000075, 0.00005, 0.000025, 0.00001, 0.0000075, 0.000005]
lapxs = np.zeros(len(dtvals))
gradxs = np.zeros(len(dtvals))
for i in range(len(dtvals)):
    laplacian = calc_laplacian(r1, r2, R, dtvals[i])
    lapxs[i] = laplacian
    grad2 = calc_grad2(r1, r2, R, dtvals[i])
    gradxs[i] = grad2[2]
    
print gradxs
    
plt.plot(-np.log(dtvals), lapxs, 'bd')
final_val = [lapxs[len(dtvals)-1]]*len(dtvals)
plt.plot(-np.log(dtvals), final_val, 'k-')
plt.ylim([-0.1886, -0.1881])
plt.yticks([-0.1886, -0.1885, -0.1884, -0.1883, -0.1882, -0.1881])
plt.xlabel('-log h')
plt.ylabel('Laplacian by central differences')
plt.show()




