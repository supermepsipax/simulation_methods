# python 3
import math
import matplotlib.pyplot as plt
import random as rnd

# Generate two Gaussian random numbers with sigma=1
def gauss_rnd(sigma):
    w = 2
    while (w >= 1):
        rx1 = 2 * rnd.random() - 1
        rx2 = 2 * rnd.random() - 1
        w = rx1*rx1 + rx2*rx2
    
    w = math.sqrt(-1*math.log(w)/w);
    return sigma*rx1*w, sigma*rx2*w

# The potential
def V(r):
    return 4*epsilon*((math.pow(sig,12)*math.pow(r,-12)) - (math.pow(sig,6)*math.pow(r,-6)))

# The force              
def F(r):
   return 4*epsilon*((12*math.pow(sig,12)*math.pow(r,-13)) - (6*math.pow(sig,6)*math.pow(r,-7)))

# Calculate the shortest periodic distance, unit cell [0,bx],[0,by]
def pbc_dist(x1,y1,x2,y2,bx,by):
    dx = x1 - x2
    dy = y1 - y2
    if dx < -0.5*bx:
        dx += bx
    if dx > 0.5*bx:
        dx -= bx
    if dy < -0.5*by:
        dy += by
    if dy > 0.5*by:
        dy -= by
    return dx, dy, math.sqrt(dx*dx + dy*dy)

# number of particles TO FILL
n = 40

# box dimension
bx = 4
by = 4

# temperature TO FILL
Temp = 300

# mass TO FILL
mass = 39.95

# time step TO FILL
dt = 0.01 
# number of time step TO FILL
nsteps =  4000

#Lennard-Jones potential constants TO FILL
sig = 0.34
epsilon=0.979
 
# Start condition
x = []
y = []
vx = []
vy = []
fx = []
fy = []
for i in range(0,n):
    # Generate initial positions separated by at least 0.35
    dmin = 0
    while dmin < 0.35:
        dmin = bx
        rx = bx*rnd.random()
        ry = by*rnd.random()
        for j in range (0,i):
            dx, dy, r = pbc_dist(x[j],y[j],rx,ry,bx,by)
            if r < dmin:
                dmin = r
            
    x.append(rx)
    y.append(ry)
    # Generate velocities from a Gaussian distribution
    KB=0.008
    factor = math.sqrt(Temp * KB / mass) 
    g1, g2 = gauss_rnd(factor)
 #   g1 = 0
 #   g2 = 0
    vx.append(g1)
    vy.append(g2)
    # Make the force array
    fx.append(0)
    fy.append(0)


outt = []
ekin = []
epot = []
etot = []
for step in range(0,nsteps+1):

    v = 0
    ek = 0

    for i in range(0,n):
        fx[i] = 0
        fy[i] = 0

    for i in range(0,n):
        for j in range(i+1,n):
            dx, dy, r = pbc_dist(x[i],y[i],x[j],y[j],bx,by)
            v  += V(r)
            fij = F(r)
            fx[i] += fij*dx/r
            fy[i] += fij*dy/r
            fx[j] -= fij*dx/r
            fy[j] -= fij*dy/r


#    ek = 0
    for i in range(0,n):
        if step > 0:
        # Update the velocities with a half step
            vx[i] += fx[i]*0.5*dt/mass
            vy[i] += fy[i]*0.5*dt/mass

        ek += 0.5*mass*(vx[i]*vx[i] + vy[i]*vy[i])

        # Update the velocities with a half step
        vx[i] += fx[i]*0.5*dt/mass
        vy[i] += fy[i]*0.5*dt/mass

        # Update the coordinates
        x[i]  += vx[i]*dt
        y[i]  += vy[i]*dt
        if x[i] < 0:
           x[i] += bx
        if x[i] >= bx:
            x[i] -= bx
        if y[i] < 0:
            y[i] += by
        if y[i] >= by:
            y[i] -= by
        
    if step % 10 == 0:
#        print step
        outt.append(step*dt)
        ekin.append(ek)
        epot.append(v)
        etot.append(v + ek)

    if step % 50 == 0:
        plt.clf()
        plt.axis([0,bx,0,by])
        plt.plot(x,y, 'ro', markersize=25)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.draw()
        plt.pause(0.02)



plt.ioff()
plt.clf()
plt.plot (outt, ekin, label = "kinetic energy")
plt.plot (outt, epot, label = "potential energy")
plt.plot (outt, etot, label = "total energy")
plt.xlabel("time $[ps]$")
plt.ylabel("Energy $[kj/mol]$")
plt.legend(loc='best', borderpad=0.5, fontsize = 10)
plt.draw()
plt.show()

