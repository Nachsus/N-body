'''
Cython time:
'''

import numpy as np 
import matplotlib.pyplot as plt  
import random
import matplotlib.animation as animation 
from multiprocessing import Process
import multiprocessing
from threading import Thread
import time

cdef double G = 1E-8
cdef double t = 4*10000 #4*10000
cdef int tsteps = 1 #4*200
cdef double tstep = t/tsteps


cdef class Body:
    cdef public double xpos
    cdef public double ypos
    cdef public double xvel
    cdef public double yvel
    cdef public double mass
    cdef public double accx
    cdef public double accy
    cdef public double xCOM 
    cdef public double yCOM 
    cdef public double afstx 
    cdef public double afsty 
    cdef public double d 
    cdef public double soft


    def __init__(self,xpos,ypos,xvel,yvel,mass):
        self.xpos = xpos
        self.ypos = ypos
        self.xvel = xvel
        self.yvel = yvel
        self.mass = mass
        self.accx = 0
        self.accy = 0


    cpdef update_vel(self):
        self.xvel += tstep * self.accx 
        self.yvel += tstep * self.accy
    

    cpdef update_pos(self):
        self.xpos += tstep * self.xvel 
        self.ypos += tstep * self.yvel
    

    cpdef add_force(self,branch):
        if branch.N > 0:
            xCOM = branch.xCOM
            yCOM = branch.yCOM
            afstx = self.xpos - xCOM
            afsty = self.ypos - yCOM
            d = np.sqrt(((afstx)**2 + (afsty)**2))
            if d != 0:
                soft = 1
                if d < soft:
                    d += soft
                    self.accx += -(G * branch.mass)/(d**2)  *  afstx/d 
                    self.accy += -(G * branch.mass)/(d**2)  *  afsty/d
                else:
                    self.accx += -(G * branch.mass)/(d**2)  *  afstx/d 
                    self.accy += -(G * branch.mass)/(d**2)  *  afsty/d


cdef class Branch:
    cdef public double xmin
    cdef public double xmax
    cdef public double ymin
    cdef public double ymax
    cdef list bodies_to_check
    cdef list bodies
    cdef public int N
    cdef list COM
    cdef double newx
    cdef double newy
    cdef object son1
    cdef object son2
    cdef object son3
    cdef object son4
    cdef list sons
    cdef public double mass
    cdef public double x
    cdef public double y
    cdef public double Mass
    cdef public double xCOM
    cdef public double yCOM
    cdef object body
    cdef public double bx 
    cdef public double by
    cdef public double l 
    cdef public double d
    

    def __init__(self,xmin,xmax,ymin,ymax,bodies_to_check):
        self.N = 0
        self.xmin = xmin
        self.xmax = xmax 
        self.ymin = ymin 
        self.ymax = ymax 
        self.bodies_to_check = bodies_to_check
        self.bodies = []
        self.bodies_in_corners()
        if self.N > 1:
            self.makeBranches()
        self.xCOM = 0
        self.yCOM = 0
        self.COM_func()
    

    cpdef between(self, bodyx, bodyy):
        if bodyx > self.xmin and bodyx < self.xmax and bodyy > self.ymin and bodyy < self.ymax:
            return True
        else:
            return False
        


    cpdef bodies_in_corners(self):
        cdef list bodies_in_area = []
        for i in self.bodies_to_check:
            if self.between(i.xpos,i.ypos):
                bodies_in_area.append(i)
        self.bodies = bodies_in_area
        self.N = len(bodies_in_area)
    

    cpdef makeBranches(self):
        newx = (self.xmin+self.xmax)/2
        newy = (self.ymin+self.ymax)/2
        son1 = Branch(self.xmin,newx,self.ymin,newy,self.bodies)
        son2 = Branch(self.xmin,newx,newy,self.ymax,self.bodies)
        son3 = Branch(newx,self.xmax,self.ymin,newy,self.bodies)
        son4 = Branch(newx,self.xmax,newy,self.ymax,self.bodies)
        sons = [son1,son2,son3,son4]
        self.sons = sons
    

    cpdef COM_func(self):
        if self.N ==0:
            self.mass = 0
        else:
            x = 0
            y = 0
            Mass = 0
            for body in self.bodies:
                x += body.xpos*body.mass
                y += body.ypos*body.mass
                Mass += body.mass 
            self.xCOM = x/Mass 
            self.yCOM = y/Mass
            self.mass = Mass
    

    cpdef theta(self,body):
        bx = body.xpos
        by = body.ypos 
        if self.xCOM != 0:
            l = np.abs(self.xmin-self.xmax)
            d = np.sqrt( (self.xCOM-bx)*(self.xCOM-bx) + (self.yCOM-by)*(self.yCOM-by) )
            if d== 0:
                return 0
            else:
                return l/d 
        else:
            return 0
    

    cpdef open(self,body):
        self.sons[0].walk(body)
        self.sons[1].walk(body)
        self.sons[2].walk(body)
        self.sons[3].walk(body)
    

    cpdef walk(self,body):
        if self.theta(body) > 0.3 and self.N > 1:
            self.open(body)
        else:
            body.add_force(self)
    



all_bodies = []
all_bodies_brute = []
N = 100000

for j in range(1):
    for i in range(int(N)):
        x = random.uniform(-90,90)
        y = random.uniform(-90,90)
        #r = 20*j
        #x = r*np.cos(2*np.pi/(N/5) * i)
        #y = r*np.sin(2*np.pi/(N/5) * i)
        m = 1
        vx = 0
        vy = 0
        '''
        if j % 2 == 0:
            vx = y/5000
            vy = -x/5000
        else:
            vx = -y/5000
            vy = x/5000
        '''
        all_bodies.append(Body(x,y,vx,vy,m))
        all_bodies_brute.append(Body(x,y,vx,vy,m))


for timestep in range(tsteps):
    start = time.time()
    br = Branch(-100,100,-100,100,all_bodies)
    for body in all_bodies:
        body.accx = 0
        body.accy = 0
        br.walk(body)
        body.update_vel()
        body.update_pos()
    slut = time.time()

print(slut-start)

