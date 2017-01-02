import math
import numpy as np

R = 4.5 #You need to match these numbers with the values you put in SC.h when you ran the data
H = 1.4

pid = [] #Where the particle id will be stored (basically useless)
sigma = [] #All the sigma values will be the same
l = [] #same with particle length
x = []
y = []
z = []
ux = []
uy = []
uz = []

INF = 'NEWFILE2' #This is the file that only contains the packing at the critical packing fraction
pid, sigma, l, x, y, z, ux, uy, uz = np.loadtxt(INF, unpack = True)

r = [] #radial distance of particle from center
dr = 0.05 #bin size
N = int(1/dr) #number of bins
n = [] #list that stores how many particles are in each bin
zavg = [] #average uz value for particles in each bin, to check if particles align with walls

def RadDist(): #turns x & y coordinates into normalized polar r coordinate
    for i in range(0,len(x)):
        r.append(math.sqrt(x[i]*x[i]+y[i]*y[i])/R)

def Sort(): #sorts the particles into their bins
    for i in range(0,N):
        a = int(0)
        b = 0
        for j in range(0,len(x)):
            if ((r[j]>=float(i)*dr) and (r[j]<=float(i+1)*dr)):
                a+=int(1)
                b+=math.sqrt(uz[j]*uz[j])
        n.append(a)
        zavg.append(b/a)
        q = (float(i)+0.5)*dr #this value is the center of the bin
        print (q, n[i], zavg[i])

def main():
    RadDist()
    Sort()

main()


