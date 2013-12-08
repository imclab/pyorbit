import numpy as np
from scipy.integrate import ode, odeint
from itertools import combinations
import pylab

#steps:
# start time -> nbody sim -> positions/velocities of main bodies
# main body sims + particle -> particle motion
# main + particle + impulse vector -> craft motion (rk4 component)    


def f(t,Y, M=None, P=None, g=6.67e-11):
    LY = len(Y)
    HY = len(Y)/2
    nbody =LY/4
    dY = np.zeros(LY)
    if not hasattr(M, '__iter__'):
        M = [1 for i in xrange(nbody)]
        
    pos = np.array(Y[:HY]).reshape((nbody,2))
    vel = np.array(Y[HY:]).reshape((nbody,2))

    dY[:HY] = np.array(Y[HY:])
    nm = np.sqrt( dY[4]**2 +  dY[5]**2) 
    
    
    if P:
        dY[4] += P* dY[4]/nm
        dY[5] += P* dY[5]/nm

    dists = np.zeros((nbody, nbody))
    diffs = np.zeros((nbody, nbody, 2))
    for c in combinations(xrange(nbody), 2):

        dists[c[0], c[1]] = sum((pos[c[0]] - pos[c[1]])**2)**1.5
        dists[c[1], c[0]] = dists[c[0], c[1]]
        
        diffs[c[0], c[1], :] = (pos[c[0]] - pos[c[1]])
        diffs[c[1], c[0], :] = -diffs[c[0], c[1], :]
    
    for i in xrange(nbody):
        xx = sum([M[k]*diffs[k,i,0]/dists[k,i] for k in xrange(nbody) if k != i])
        yy = sum([M[k]*diffs[k,i,1]/dists[k,i] for k in xrange(nbody) if k != i])
        dY[HY + 2*i] = g * xx
        dY[HY + 2*i + 1] = g * yy

    return dY

# all units: m, k, s

M = [5.97219e24, 7.3477e22 ,10]
pos = [[0,0],
       [405503000,0],
        [-35786000,0]]
vel = [[0, 0],
       [0.,964.],
        [0,-3430]]

dt = 60*60*1
t0 = 0
tfinal = 60*60*24*27 # final start time
nbody = len(M)
y0 = [val for subl in pos for val in subl] + [val for subl in vel for val in subl]

T = np.arange(0,tfinal,dt)

#solution = odeint(f, y0,T,args=(M,P))

solution = ode(f).set_integrator('dopri5',atol=1e-6)
P=[0,0]
solution.set_initial_value(y0,t0).set_f_params(M, P)

cols = ['b','r','g']


pylab.figure()
for x,y in pos:
    pylab.plot(x,y,'+',markersize=10)
    
k = 0
while solution.successful() and solution.t < tfinal: 
    if k == 90:    
        solution.set_f_params(M, 15000)
    else:
        solution.set_f_params(M, 0)
    y = solution.integrate(solution.t+dt)
    for i in xrange(nbody):
        pylab.plot(y[2*i], y[2*i+1], cols[i]+'.',markersize=2)
    k+=1
pylab.show() # show plot