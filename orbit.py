import numpy as np
from scipy.integrate import ode, odeint
from itertools import combinations
import pylab


class Orbits(object):

    """
    Class models the gravitational dynamics of n-bodies,
    along with a specified number of small craft (ships, satellites, etc.)

    Craft have negligible mass compared to the bodies, and can have impulse
    applied at any moment of time to change orbits or move between bodies.
    """

    def __init__(self):
        self.g = 6.67e-11
        self.M = []
        self.T = []
        self.t0 = 0
        self.dt = 60 * 60
        self.tfinal = 60 * 60 * 24 * 28
        self.T = np.arange(0, self.tfinal, self.dt)
        self.ntimes = len(self.T)
        self.pos = []
        self.vel = []
        self.nbody = 0

        self.ncraft = 0
        self.cvel = []
        self.cpos = []
        
        self.impulses=[]

        self.body_solutions = []

    def set_dt(self, dt):
        """
        Sets timestep for ODEs
        """
        self.dt = dt
        self.T = np.arange(0, self.tfinal, dt)
        self.ntimes = len(self.T)

    def set_tfinal(self, t):
        """
        Sets final time for ODEs
        """
        self.tfinal = t
        self.T = np.arange(0, tfinal, dt)
        self.ntimes = len(self.T)

    def add_body(self, m, pos0, vel0):
        """
        Adds body (planet, star, moon, etc) to problem.
        Must specify mass, initial position, and velocity.
        """
        self.nbody += 1
        self.M.append(m)
        self.pos.append(pos0)
        self.vel.append(vel0)

    def add_craft(self, pos0, vel0):
        """
        Adds craft (object with mass negligible compared to bodies)
        to problem. must specify initial position and velocity.
        """
        self.ncraft += 1
        self.cpos.append(pos0)
        self.cvel.append(vel0)
        dv = np.zeros(self.ntimes)
        self.impulses.append(dv)

    def solve_bodies(self):
        """
        Solves body positions+velocities
        """
        y0 = [val for subl in self.pos for val in subl] + \
            [val for subl in self.vel for val in subl]

        #solution = ode(self.f).set_integrator('dopri5', atol=1e-6)
        #solution.set_initial_value(y0, self.t0).set_f_params(self.M)

        self.body_solutions = odeint(self.f, y0, self.T)

    def solve_crafts(self):
        """
        Solves craft positions+velocities
        """
        self.craft_solutions = np.zeros((self.ncraft, self.ntimes, 4))

        for j in xrange(self.ncraft):
            if not len(self.body_solutions):
                self.solve_bodies()
            y0 = self.cpos[j] + self.cvel[j]
            self.craft_solutions[j,0] = y0
            solution = ode(self.f_craft).set_integrator('lsoda', rtol=1e-6)
            solution.set_initial_value(y0, self.t0)#.set_f_params(self.M)
            i = 1
            for t in self.T[1:]:
                solution.set_f_params(self.body_solutions[i], self.impulses[j][i])
                y = solution.integrate(t)
                self.craft_solutions[j,i] = y
                i+=1

    def plot_bodies(self):
        """
        Plot body positions
        """
        cols = ['b', 'r', 'g']

        pylab.figure()
        for x, y in self.pos:
            pylab.plot(x, y, '+', markersize=10)

        for k in xrange(self.ntimes):
            y = self.body_solutions[k]
            for i in xrange(self.nbody):
                pylab.plot(y[2 * i], y[2 * i + 1], cols[i] + '.', markersize=2)

    def plot(self):
        """
        Plot body and craft positions
        """
        bcols = ['b', 'r', 'g']
        ccols = ['g','y','k']
        pylab.figure()
        for x, y in self.pos:
            pylab.plot(x, y, '+', markersize=10)

        for k in xrange(self.ntimes):
            y = self.body_solutions[k]
            for i in xrange(self.nbody):
                pylab.plot(y[2 * i], y[2 * i + 1], bcols[i] + '.', markersize=2)
                
            for i in xrange(self.ncraft):
                y = self.craft_solutions[i,k]
                pylab.plot(y[0], y[1], ccols[i] + '.', markersize=2)

    def f_craft(self, t, Y, b_sol=None, dV=None):
        """RHS of single craft position and velocity ODE.
        Body positions and velocities must be pre-solved."""
        
        LY = len(b_sol)
        HY = LY / 2
        pos = np.array(b_sol[:HY]).reshape((self.nbody, 2))
        dY = np.zeros(4)
        diffsx = np.zeros(self.nbody)
        diffsy = np.zeros(self.nbody)
        dists = np.zeros(self.nbody)
        dY[:2] = Y[2:]
        if dV:
            #dV = dV*self.dt/1000.      
            nm = np.linalg.norm(dY[:2], 2)
            dY[0] += dV * dY[0] / nm
            dY[1] += dV * dY[1] / nm

        for i in xrange(self.nbody):
            diffsx[i] = Y[0] - pos[i][0]
            diffsy[i] = Y[1] - pos[i][1]
            dists[i] = sum((Y[:2] - pos[i]) ** 2) ** 1.5

        dY[2] = self.g*sum([-self.M[i]*diffsx[i]/dists[i] for i in xrange(self.nbody)])
        dY[3] = self.g*sum([-self.M[i]*diffsy[i]/dists[i] for i in xrange(self.nbody)])
        return dY

    def f(self, Y, t):
        """RHS of body position ODE"""

        M = self.M
        LY = len(Y)
        HY = LY / 2
        nbody = self.nbody
        dY = np.zeros(LY)

        pos = np.array(Y[:HY]).reshape((nbody, 2))

        dY[:HY] = np.array(Y[HY:])
        nm = np.sqrt(dY[4] ** 2 + dY[5] ** 2)

        dists = np.zeros((nbody, nbody))
        diffs = np.zeros((nbody, nbody, 2))
        for c in combinations(xrange(nbody), 2):

            dists[c[0], c[1]] = sum((pos[c[0]] - pos[c[1]]) ** 2) ** 1.5
            dists[c[1], c[0]] = dists[c[0], c[1]]

            diffs[c[0], c[1], :] = (pos[c[0]] - pos[c[1]])
            diffs[c[1], c[0], :] = -diffs[c[0], c[1], :]

        for i in xrange(nbody):
            xx = sum([M[k] * diffs[k, i, 0] / dists[k, i]
                     for k in xrange(nbody) if k != i])
            yy = sum([M[k] * diffs[k, i, 1] / dists[k, i]
                     for k in xrange(nbody) if k != i])
            dY[HY + 2 * i] = self.g * xx
            dY[HY + 2 * i + 1] = self.g * yy

        return dY

if __name__ == "__main__":
    # Set up system, will be Earth-centered inertial frame (ECI)
    sys = Orbits()
    sys.dt = 60*60
    # Add the Earth and Moon (m-k-s units)
    sys.add_body(5.97219e24, [0, 0], [0, 0])  # Earth
    sys.add_body(7.3477e22, [405503000, 0], [0, 964.])  # Moon
    
    # Place a craft in a geostationary orbit (~35,000 km) 
    sys.add_craft([-35786000, 0], [0, -3340])  
    # Place small craft into a lower orbit
    sys.add_craft([-4140000, 0], [0, -7700]) 
    
    # apply some impulse (kind of unstable at the moment...)
    # boost first sattellite towards moon
    sys.impulses[0][350] = 16000 
    # boost 2nd to geostationary
    sys.impulses[1][75] = 3500 
    
    sys.solve_crafts()
    
    
    sys.plot()
    pylab.show()
