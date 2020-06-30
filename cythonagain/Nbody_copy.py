import numpy as np 
import matplotlib.pyplot as plt  
import random
import matplotlib.animation as animation 
import time

from multiprocessing import Process
import multiprocessing

from threading import Thread


#G = 6.67*10**(-11)
G = 10**(-8)
t = 4*10000 #4*10000
tsteps = 1 #4*200
tstep = t/tsteps


'''
The body class is just the bodies, with some position and velocity.
'''
class body:
    def __init__(self,pos,vel,mass):
        self.pos = pos 
        self.vel = vel
        self.mass = mass
        self.accx = 0
        self.accy = 0
        self.accxbrute = 0
        self.accybrute = 0
        self.counter = 1
    

    def update_vel(self):
        dvx = tstep * self.accx
        dvy = tstep * self.accy
        self.vel = [self.vel[0]+dvx, self.vel[1]+dvy]
    

    def update_pos(self):
        dx = tstep * self.vel[0]
        dy = tstep * self.vel[1]
        self.pos = [self.pos[0]+dx, self.pos[1]+dy]
    
    
    def add_force(self,branch):
        if branch.N > 0:
            xCOM = branch.COM[0]
            yCOM = branch.COM[1]
            afstx = self.pos[0] - xCOM
            afsty = self.pos[1] - yCOM
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
    

    def add_force_brute(self,body2): #for testig purposes
        afstx = self.pos[0] - body2.pos[0]
        afsty = self.pos[1] - body2.pos[1]
        d = np.sqrt((afstx)**2 + (afsty)**2)
        if d == 0:
            print('hello')
        if d < 1:
            d += 1
            self.accxbrute += -(G * body2.mass)/(d**2)  *  afstx/d
            self.accybrute += -(G * body2.mass)/(d**2)  *  afsty/d
            self.counter += 1
        else:
            self.accxbrute += -(G * body2.mass)/(d**2)  *  afstx/d
            self.accybrute += -(G * body2.mass)/(d**2)  *  afsty/d
            self.counter += 1


    


# Making some bodies just to see if the code works
all_bodies = []
all_bodies_brute = []
N = 1000

for j in range(1):
    for i in range(int(N/5)):
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
        all_bodies.append(body([x,y],[vx,vy],m))
        all_bodies_brute.append(body([x,y],[vx,vy],m))


'''
This is the branch class. It makes up boxes in the system which contains bodies.
If this branch contains more than one body, it makes 4 new branches. 
It does this untill there is maximum 1 body in each box.
'''
class branch:
    def __init__(self, corners, bodies_to_check, father, gen):
        self.corners = corners
        self.N = 0 # Number of bodies in branch
        self.sons = [] # Son branches
        self.bodies = []
        self.bodies_to_check = bodies_to_check
        self.father = father
        self.gen = gen
        self.bodies_in_corners()
        if self.N > 1:
            self.makebranches()
        
        self.COM = []
        self.COM_func()
        

    def between(self, bodyposition): # Simply finds out if the bodyposition is within the limits
        if bodyposition[0] > self.corners[0][0] and bodyposition[0] < self.corners[0][1] and bodyposition[1] > self.corners[1][0] and bodyposition[1] < self.corners[1][1]:
            return True
        else:
            return False


    def bodies_in_corners(self): # Makes a list of all bodies inside the limits
        bodies_in_area = []
        for i in range(len(self.bodies_to_check)):
            if self.between(self.bodies_to_check[i].pos):
                bodies_in_area.append(self.bodies_to_check[i])
        self.bodies = bodies_in_area
        self.N = len(bodies_in_area)

    
    def makebranches(self): # Makes 4 new branches inside the old one with new limits.
        newx = (self.corners[0][0] + self.corners[0][1]) / 2
        newy = (self.corners[1][0] + self.corners[1][1]) / 2
        xmin = self.corners[0][0]
        xmax = self.corners[0][1]
        ymin = self.corners[1][0]
        ymax = self.corners[1][1]
        gen = (self.gen+1)
        son1 = branch([[xmin,newx],[ymin,newy]],self.bodies,self,gen)
        son2 = branch([[xmin,newx],[newy,ymax]],self.bodies,self,gen)
        son3 = branch([[newx,xmax],[ymin,newy]],self.bodies,self,gen)
        son4 = branch([[newx,xmax],[newy,ymax]],self.bodies,self,gen)
        sons = [son1,son2,son3,son4]
        self.sons = sons
    

    def COM_func(self): # Finds the COM of the box
        if self.N == 0:
            self.mass = 0
        else:
            x = 0
            y = 0
            Mass = 0
            for i in self.bodies:
                x += (i.pos[0]*i.mass)
                y += (i.pos[1]*i.mass)
                Mass += i.mass
            COMx = x/Mass
            COMy = y/Mass
            self.mass = Mass 
            self.COM = [COMx, COMy]
    

    def theta(self,body):
        bx = body.pos[0]
        by = body.pos[1]
        if len(self.COM) != 0:
            xCOM = self.COM[0]
            yCOM = self.COM[1]
            l = np.abs(self.corners[0][0]-self.corners[0][1])
            d = np.sqrt(((xCOM-bx)**2 + (yCOM-by)**2))
            if d == 0:
                return 0
            else:
                return l/d
        else:
            return 0
    

    def open(self, body):
        self.sons[0].walk(body)
        self.sons[1].walk(body)
        self.sons[2].walk(body)
        self.sons[3].walk(body)


    def walk(self,body):
        if self.theta(body) > 0.3 and self.N > 1:
            self.open(body)
        else:
            body.add_force(self)


firstcorns = [[-100,100],[-100,100]]

'''
arr = [2345,547,679,4567,2345,2345,457,698,709,5679,4567,3456,2345,437,5679,7608,879,789]
arr_arr = []
for i in range(6):
    arr_arr.append(arr[i:i+int(len(arr)/6)])

print(arr_arr)
'''
'''
body opdeling test
'''
manager = multiprocessing.Manager()

body_lists = manager.list()
cores = 1
for i in range(cores):
    body_lists.append(all_bodies[i:i + int(len(all_bodies)/cores)])

#print(len(body_lists[0]))
dembodies = manager.list()
def sims(dembodies):
    arr = []
    for body in dembodies:
        body.accx = 0
        body.accy = 0
        Node = br
        Node.walk(body)
        body.update_vel()
        body.update_pos()
        #pos.append(body.pos)
        arr.append(body)
    dembodies = arr


#ttree1 = time.time()
#print(body_lists[1])


thred1 = time.time()
print(all_bodies[0].pos)
print(all_bodies[-1].pos)
body_pos_t = []
for time_step in range(tsteps):
    tid1 = time.time()
    #manager = multiprocessing.Manager()
    #pos = manager.list()
    #print(body_lists[0][0].pos)
    br = branch(firstcorns,all_bodies,None,1)
    for body in all_bodies:
        body.accx = 0
        body.accy = 0
        Node = br
        Node.walk(body)
        body.update_vel()
        body.update_pos()
        #pos.append(body.pos)
        #arr.append(body)
    #p0 = Process(target=sims, args=([all_bodies]))
    #p1 = Process(target=sims, args=([body_lists[1]]))
    #p2 = Process(target=sims, args=([body_lists[2]]))
    #p3 = Process(target=sims, args=([body_lists[3]]))
    #p4 = Process(target=sims, args=([body_lists[4]]))
    #p5 = Process(target=sims, args=([body_lists[5]]))
    #p6 = Process(target=sims, args=([body_lists[6]]))
    #p7 = Process(target=sims, args=([body_lists[7]]))


    #p0.start()
    #p1.start()
    #p2.start()
    #p3.start()
    #p4.start()
    #p5.start()
    #p6.start()
    #p7.start()

    #p0.join()
    #p1.join()
    #p2.join()
    #p3.join()
    #p4.join()
    #p5.join()
    #p6.join()
    #p7.join()
    
    #print(pos)
    #print(poses)
    #for body in body_lists[0]:
    #    pos.append(body.pos)

    #print('--'*10)
    #print(pos)
    #print('--'*10)
    #body_pos_t.append(pos)
    tid2 = time.time()
    print(tid2-tid1)

thred2 = time.time()

print(all_bodies[0].pos)
print(all_bodies[-1].pos)
'''
nothred1 = time.time()
for time_step in range(tsteps):
    for body in all_bodies:
            body.accx = 0
            body.accy = 0
            Node = br
            Node.walk(body)
            body.update_vel()
            body.update_pos()
nothred2 = time.time()
'''
#print(thred2-thred1)
#print(nothred2-nothred1)



'''
body_pos_t = []
count = 1
for time_step in range(tsteps):
    #print(count/tsteps)
    #count += 1
    br = branch(firstcorns,all_bodies,None,1)
    poss = []
    for body in all_bodies:
        body.accx = 0
        body.accy = 0
        Node = br
        Node.walk(body)
        body.update_vel()
        body.update_pos()
        poss.append(body.pos)
    body_pos_t.append(poss)
'''
#ttree2 = time.time()

#print(ttree2-ttree1)

'''
tbrute1 = time.time()
for time_step in range(tsteps):
    #count = 1
    for body in all_bodies_brute:
        #count += 1
        #print(count/len(all_bodies_brute))
        body.accxbrute = 0
        body.accybrute = 0
        for body2 in all_bodies_brute:
            if body != body2:
                body.add_force_brute(body2)
tbrute2 = time.time()

errsx = []
errsy = []
for i in range(len(all_bodies)):
    errsx.append(all_bodies[i].accx/all_bodies_brute[i].accxbrute)
    errsy.append(all_bodies[i].accy/all_bodies_brute[i].accybrute)
'''



'''
plt.figure()
plt.title('Box')
count = 1
for body in all_bodies:
    print(count)
    count += 1
    plt.scatter(body.pos[0],body.pos[1], s=1)

plt.figure()
plt.title('Brute')
count = 1
for body in all_bodies_brute:
    print(count)
    count += 1
    plt.scatter(body.pos[0],body.pos[1], s=1)
'''
'''
plt.figure()
for i in range(len(all_bodies)):
    f1 = abs(all_bodies_brute[i].accxbrute + all_bodies_brute[i].accybrute)/2
    f2 = abs(all_bodies[i].accx + all_bodies[i].accy)/2
    x = f1
    y = abs(f2-f1)/f1
    plt.scatter(x,y)


print('Tider ---')
print(ttree2-ttree1)
print(tbrute2-tbrute1)
print(((ttree2-ttree1)/(tbrute2-tbrute1))**(-1))
print('----')
'''
#print(errsx)
#print(sum(errsx)/N)
#print(errsy)
#print(sum(errsy)/N)

#print((sum(errsx)+sum(errsy))/(2*N))

#plt.show()


'''
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1,1,1)
u = 1

def ani_scatter(x):
    l = x*u
    print(l/(tsteps/u)/u)
    ax.clear()
    ax.grid()
    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    for step in body_pos_t[l]:
        plt.scatter(step[0],step[1],s=1)

x_tot = []
y_tot = []
for time in body_pos_t:
    x = []
    y = []
    x.append(-100)
    y.append(-100)
    x.append(100)
    y.append(-100)
    x.append(100)
    y.append(100)
    x.append(-100)
    y.append(100)
    for pos in time:
        if pos[0] < 100 and pos[0] > -100:
            if pos[1] < 100  and pos[1] > -100:
                x.append(pos[0])
                y.append(pos[1])
    x_tot.append(x)
    y_tot.append(y)

import seaborn as sns


def ani_hist(x):
    l = x*u
    print(l/(tsteps/u)/u)
    ax.clear()
    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    h, x1, y1, p = plt.hist2d(x_tot[l],y_tot[l],bins=20)
    ax.clear
    plt.imshow(h, origin='lower', interpolation = 'gaussian')
    #sns.kdeplot(x_tot[l],y_tot[l])


anim = animation.FuncAnimation(fig,ani_scatter,frames=int(tsteps/u),interval=40)

anim.save('test.mp4')
#plt.show()
'''
'''
In general i see that my box code give an acceleration well within 1% of a normal brute force method.
Though sometimes i see a variense in around 30%. I don't know why.
'''



'''
mpi4py
cython
https://github.com/jmd-dk 
https://arxiv.org/abs/1510.07621 

def evaluate_force:
  for body in bodies:
    node = root_node
    if ok_precision(body,node):
      dd_force(body, node)
    else:
      pen_node(body,node)

def open_node(body, node):
  for nodes in node:
    if ok(body,node):
      add_force(body, node)
    else
      open_node(body, node)
def evaluate_force:
  for body in bodies:
    node = root_node
    if ok_precision(body,node):
      dd_force(body, node)
    else:
      pen_node(body,node)

def open_node(body, node):
  for nodes in node:
    if ok(body,node):
      add_force(body, node)
    else
      open_node(body, node)
'''