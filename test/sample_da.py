"""
referance:
https://qiita.com/elect-gombe/items/82b60123775c178258f4#letkf
"""


import time
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import pandas as pd
import sys
import json
import numba as nb
import Bmat
import os
import scipy


F=8
J=40
ADAY=0.2
AYEAR=365
YEARS_TO_SIMULATE=2
dt=0.05
t_spinup=ADAY*AYEAR*1
tend=AYEAR*ADAY*YEARS_TO_SIMULATE
t=np.arange(0,tend+dt,dt)
ft=np.arange(ADAY*AYEAR,tend+dt,dt)

x=np.zeros(J)
v=np.zeros(J)
RMSE=[]
xgtlist=[]
xobzlist=[]
STEP=2
FIRSTN=40
delta=0.1
H=np.identity(J)[:FIRSTN:STEP]
R=np.identity(int(FIRSTN/STEP))

printB=False

#use 5 points formula
# @nb.jit
# def differential(x,h,hr,dt):
#    return (1.0/hr)*(
#        -(1./12.)*runge4th(x+2*h,dt)
#        +(2./3.)*runge4th(x+h,dt)
#        -(2./3.)*runge4th(x-h,dt)
#        +(1./12.)*runge4th(x-2*h,dt)
#         )

@nb.jit
def differential(x,h,hr,dt):
   return (runge4th(x+h,dt)-runge4th(x,dt))/hr

@nb.jit
def invnumbajit(A):
    return np.linalg.inv(A)
@nb.jit
def solveacc(x):
    a=np.empty(J)
    for j in range(J):
        a[j]=(x[(j+1)%J]-x[(j-2+J)%J])*x[(j-1+J)%J]-x[j]+F
    return a


@nb.jit
def euler(x,dt):
    a=solveacc(x)
    return x+a*dt

@nb.jit
def runge4th(x,dt):
    a=solveacc(x)
    k=np.empty((5,J))
    k[1]=a*dt
    a=solveacc(x+k[1]*0.5)

    k[2]=a*dt
    a=solveacc(x+k[2]*0.5)

    k[3]=a*dt
    a=solveacc(x+k[3])

    k[4] = a*dt
    k[0]=(k[1]+2*k[2]+2*k[3]+k[4])/6.
    return x+k[0]


def initialize():
    global x
    x=np.ones(J)*F
    x[int(J/2)+1]*=1.001

def rmse(x):
    return np.sqrt(np.average(
        np.square(x)
    ))

def evaluate_error():
    rmsevlist=np.zeros(100)
    for j in range(10000):
        x=xgtlist[1000]+np.random.randn(J)/10000
        for i,xnow in enumerate(xgtlist[1001:1100]):
            x=runge4th(x,dt)
            diff=x-xnow
            rmsev=rmse(diff)
            rmsevlist[i]+=(rmsev)

    plt.title("RMSE")
    plt.yscale('log')
    plt.plot(t[:99],rmsevlist[:99]/rmsevlist[0])
    plt.show()

def runsimulation(spinup):
    global x
    initialize()
    for i in range(spinup):
        x=runge4th(x,dt) #Update x

    for i,tnow in enumerate(t):
        x=runge4th(x,dt) #Update x
        xgtlist.append(np.copy(x))

def loaddata(path):
    global F,J,ADAY,AYEAR,YEARS_TO_SIMULATE,dt,t_spinup,tend,xobzlist,xgtlist
    observation_path=path+"_obz.npy"
    groundtruth_path=path+"_gt.npy"
# load parameter
    with open(path+"cfg.json") as f:
        df = json.load(f)
        F,J,ADAY,AYEAR=df["F"],df["J"],df["ADAY"],df["AYEAR"]
        YEARS_TO_SIMULATE,dt,t_spinup,tend=df["YEARS_TO_SIMULATE"],df["dt"],df["t_spinup"],df["tend"]
        observation_path,groundtruth_path=df["observation_filepath"],df["groundtruth_filepath"]

# load experiment data
    xobzlist=np.load(observation_path)
    xgtlist=np.load(groundtruth_path)

def storedata(path):
    observation_path=path+"_obz.npy"
    groundtruth_path=path+"_gt.npy"
    param={"F":F,
           "J":J,
           "ADAY":ADAY,
           "AYEAR":AYEAR,
           "YEARS_TO_SIMULATE":YEARS_TO_SIMULATE,
           "dt":dt,
           "t_spinup":t_spinup,
           "tend":tend,
           "observation_filepath":observation_path,
           "groundtruth_filepath":groundtruth_path,
    }

    with open(path+"cfg.json","w") as f:
        json.dump(param,f)

    np.save(groundtruth_path,np.array(xgtlist))
    for i in xgtlist:
        xobzlist.append(i+np.random.randn(len(i)))
    np.save(observation_path,np.array(xobzlist))

def showmat(m):
   plt.imshow(m,interpolation='nearest',cmap='jet')
   plt.colorbar()
   plt.show()

def localizescale(sigma,i):
    r=abs(i)
    r=min(r,abs(J-r))
    return np.exp(-r*r/(2*sigma*sigma))


def localizematrix(sigma,H):
    m=[]
    for i in H:
        idx=(np.where(i!=0)[0][0])
        l=[]
        for j in range(J):
            l.append(localizescale(sigma,j-idx))
        m.append(np.array(l))
    return np.array(m)

def calcinvR(sigma,idx):
    nobs=H.shape[0]
    mat=np.identity(nobs)
    for idx2,j in enumerate(mat):
        idx2=(np.where(H[idx2]!=0)[0][0])
        j*=(localizescale(sigma,idx-idx2))
    return mat

def printmat(B):
    for i in B:
        s="["
        for j in i:
            s+=("%e, "%j)
        print(s+"],")


def observation_densitymatrix_turn(n):
    delnum = J-n
    dellist=np.linspace(0,J,delnum+2) #delete grid point list. delete except first and last point.
    matt = list(np.identity(J))
    for idx,i in enumerate(dellist[1:-1]):
        matt.pop(int(i-idx))
    return np.array(matt)

def observation_densitymatrix_dense(n):
    return np.identity(J)[:n]

plotlist=[[],[],[]]
def threedvar(scale,odm):
    B=np.array(Bmat.B)*scale
    xa=xobzlist[0]
    xf=xgtlist
    IMAT=invnumbajit(np.dot(np.dot(H,B),H.T)+R)
    K=np.dot(np.dot(B,H.T),IMAT)
    for idx,i in enumerate(ft):
        xf=runge4th(xa,dt)

        xa = xf + np.dot(K,np.dot(odm,xobzlist[idx])-np.dot(H,xf))
        plotlist[0].append(rmse(xa-xgtlist[idx]))

    return np.mean(plotlist[0][500:])


def KalmanFilter(delta,odm):
    Pflist=[]
    Pa = np.identity(J)*5
#    print(xgtlist[np.random.randint(int(ADAY*400/dt),int(tend/dt))])
    xa=xobzlist[0]
    xf=xgtlist

    JM=np.empty((J,J))
    for idx,i in enumerate(ft):
        xf=runge4th(xa,dt)

        for idx2,j in enumerate(np.identity(J)):
            JM[idx2]=differential(xa,j*0.001,0.001,dt)
        Pf=(1+delta)*np.dot(np.dot(JM.T,Pa),JM)

        K=np.dot(np.dot(Pf,H.T),invnumbajit(np.dot(np.dot(H,Pf),H.T)+R))
        xa = xf + np.dot(K,np.dot(odm,xobzlist[idx])-np.dot(H,xf))
        Pa = np.dot((np.identity(J)-np.dot(K,H)),Pf)
        Pflist.append(np.copy(Pf))

        plotlist[0].append(rmse(xa-xgtlist[idx]))
        plotlist[1].append(rmse(xgtlist[idx]-xobzlist[idx]))
        plotlist[2].append(np.sqrt(np.trace(Pa)/J))

    average=Pflist[0]
    for i in Pflist[-200:]:
        average+=i
    average /= 200

    if printB==True:
        printmat(average)

    return np.mean(plotlist[0][500:])


def LETKF(delta,odm,nmem=12):
    mask=localizematrix(3,H)
    Pa = np.identity(J)*5
#    print(xgtlist[np.random.randint(int(ADAY*400/dt),int(tend/dt))])
    xa=np.array([xobzlist[0] for i in range(nmem)])+np.random.normal(0,1,(nmem,J))
    xanext=np.empty(nmem)
    plotlist=[[],[],[]]
    Palst=[]
    JM=np.empty((J,J))
    sdelta=np.sqrt(1+delta)
    for idx,i in enumerate(ft):
        xf=np.array([runge4th(i,dt) for i in xa])
        xfa=np.mean(xf,axis=0)
        dxf=np.array([i-xfa for i in xf])*sdelta
        yfa=np.mean(np.dot(H,xf.T),axis=1)
        dyf=np.array([np.dot(H,i)-yfa for i in xf])*sdelta
        xanext=[]
        for j in range(J):
            invR=calcinvR(3,j).T
            C=np.dot(dyf,invR)
            w,v=np.linalg.eig(np.identity(nmem)*(nmem-1)+np.dot(C,dyf.T))
            w=np.real(w)
            v=np.real(v)
            p_invsq=numpy.diag(1/np.sqrt(w))
            p_inv=numpy.diag(1/w)
            Wa=v @ p_invsq @ v.T
            Was=v @ p_inv @ v.T

            xanext.append((np.matlib.repmat(xfa,nmem,1)+np.dot(dxf.T,
                     np.linalg.multi_dot([Was,C,
                                          np.linalg.multi_dot([H,np.matlib.repmat(xobzlist[idx]-xfa,nmem,1).T])
                     ])+np.sqrt(nmem-1)*Wa).T).T[j])
        xa=np.array(xanext).T
        rmsev=rmse(np.mean(xa,axis=0)-xgtlist[idx])
        if idx&0x1F==0:
            sys.stdout.write(".")
            sys.stdout.flush()
        plotlist[0].append(rmsev)
        Pa=np.dot(dxf.T,dxf)/(nmem-1)
        Palst.append(Pa)
        plotlist[1].append(np.sqrt(np.trace(Pa)/J))

    average=Palst[0]
    for i in Palst[-200:]:
        average+=i
    average /= 200

#    showmat(mask*average)

    # plt.plot(plotlist[0],label="rmse")
    # plt.plot(plotlist[1],label="trace pa")
    # plt.legend()
    # plt.show()
    return np.mean(plotlist[0][500:])

def PO(delta,odm,nmem=16):
    mask=localizematrix(3,H)
    Pa = np.identity(J)*5
#    print(xgtlist[np.random.randint(int(ADAY*400/dt),int(tend/dt))])
    xa=np.array([xobzlist[0] for i in range(nmem)])+np.random.normal(0,1,(nmem,J))
    plotlist=[[],[],[]]
    Palst=[]
    JM=np.empty((J,J))
    sdelta=1#delta
    for idx,i in enumerate(ft):
        xf=np.array([runge4th(i,dt) for i in xa])
        xfa=np.mean(xf,axis=0)
        dxf=np.array([i-xfa for i in xf])*sdelta
        yfa=np.mean(np.dot(H,xf.T),axis=1)
        dyf=np.array([np.dot(H,i)-yfa for i in xf])*sdelta
        K=mask.T*np.linalg.multi_dot([dxf.T,invnumbajit(np.identity(nmem)+np.linalg.multi_dot([dyf,R,dyf.T])),dyf,R])/np.sqrt(nmem-1)
        t=np.dot(H,np.array(np.matlib.repmat(xobzlist[idx],nmem,1)+np.random.normal(0,1,(nmem,J))*(1+delta)-xf).T)
        xa=xf+np.dot(K,t).T
        plotlist[0].append(rmse(np.mean(xa,axis=0)-xgtlist[idx]))
        Pa=np.dot(dxf.T,dxf)/(nmem-1)
        Palst.append(Pa)
        plotlist[1].append(np.sqrt(np.trace(Pa)/J))

    average=Palst[0]
    for i in Palst[-200:]:
        average+=i
    average /= 200

    # showmat(mask*average)

    # plt.plot(plotlist[0],label="rmse")
    # plt.plot(plotlist[1],label="trace pa")
    # plt.legend()
    # plt.show()
    return np.mean(plotlist[0][500:])

def main():
    global x,H,R,printB
    plt.grid(which = "major", axis = "x", color = "blue", alpha = 0.8, linestyle = "--", linewidth = 1)
    plt.grid(which = "major", axis = "y", color = "green", alpha = 0.8, linestyle = "--", linewidth = 1)
    plt.xlabel('observation data count')
    plt.ylabel('RMSE')

# if you want simulation data, uncomment and run simulation.
    # print("running simulation"+str(time.time()-begintime))
    # runsimulation(int(AYEAR*ADAY/dt))
    # print("storing data"+str(time.time()-begintime))
    # storedata("data/result")
    loaddata("data/result")

    rmselst=[]
    nlist=np.arange(10,40+1,1)
    deltaset=np.arange(0.04,0.4,0.01)
    deltaset2=np.arange(0.,0.2,0.01)
    scaleset=np.arange(0.1,2,0.2)
    obzspace=(observation_densitymatrix_turn,observation_densitymatrix_dense)
    filters=(LETKF,PO,KalmanFilter,threedvar)
    for densmat in obzspace:
        print("dens:%s"%densmat.__name__)
        for fil in filters:
            print("filters:%s"%fil.__name__)
            if fil == threedvar:
                params=scaleset
            elif fil == LETKF:
                params=deltaset2
            else:
                params=deltaset

            rmseminlst=[]
            for n in nlist:
                rmselst=[]
                print("obzdata :%d"%n)
                matt=densmat(n)
                R=np.identity(n)
                H=matt
                #try finding minimum rmse for each deltaset
                for i in params:
                    if True: # any error to be caught and fail if you set False(debug only).
                        try:
                            rmse=fil(i,matt)
                        except:
                            # 999 shows error
                            rmse=999
                    else:
                        rmse=fil(i,matt)
                    pair=(i,rmse)
                    rmselst.append(rmse)
                    print("param:%3f,rmse:%3f"%pair)
                    plotlist[0],plotlist[1],plotlist[2]=[],[],[]

                print("finished in "+str(time.time()-begintime))
                rmseminlst.append(min(rmselst))

            plt.plot(nlist,rmseminlst,label="%s"%fil.__name__)

        plt.grid(which = "major", axis = "x", color = "blue", alpha = 0.8, linestyle = "--", linewidth = 1)
        plt.grid(which = "major", axis = "y", color = "green", alpha = 0.8, linestyle = "--", linewidth = 1)
        plt.xlabel('observation data count')
        plt.ylabel('RMSE')
        plt.title("minimum RMSE")
        plt.ylim(0,2)
        plt.legend()
        plt.show()


if __name__ == '__main__':
    main()