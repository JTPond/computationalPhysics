#! /usr/bin/python
import numpy as np
import itertools
from matplotlib import pyplot as plt

def get_states(N):
    out = []
    for i in itertools.combinations_with_replacement([0,1],N):
        for j in itertools.permutations(i,N):
            o = list(j)
            if not o in out:
                out.append(o)
    return out #returns list of all possible states of N spins

def get_diag_H(state):
    ran = range(len(state))
    out = 0.0
    for i in ran:
        if i != ran[-1]:
            if state[i] == state[i+1]:
                out += -1.0
            else:
                out += 1.0
        else:
            if state[i] == state[0]:
                out += -1.0
            else:
                out += 1.0
    return out

def get_off_diag_H(state1,state2,hoJ):
    flips = 0
    out = 0.0
    for i in range(len(state1)):
        if state1[i] == state2[i]:
            flips +=1
    if flips > 1 or flips < 1:
        return 0.0
    else:
        return hoJ

def get_H_hat(N,hoJ):
    ret = np.zeros(shape=(2**N,2**N))
    states = get_states(N)
    for i in range(2**N):
        for j in range(2**N):
            if i == j:
                ret[i,j] = get_diag_H(states[i])
            else:
                ret[i,j] = get_off_diag_H(states[i],states[j],hoJ)
    return ret

def E(e_vals,T,N):
    topTerm = 0.0
    botTerm = 0.0
    for em in e_vals:
        topTerm += em*np.e**(-em/T)
        botTerm += np.e**(-em/T)
    return (1.0/float(N))*(topTerm/botTerm)

def C(e_vals,T,N):
    topTerm1 = 0.0
    topTerm2 = 0.0
    botTerm = 0.0
    for em in e_vals:
        topTerm1 += em*np.e**(-em/T)
        topTerm2 += em*em*np.e**(-em/T)
        botTerm += np.e**(-em/T)
    return (1.0/N)*(1.0/(T*T))*(-((topTerm1**2.0)/(botTerm**2.0)) + (topTerm2/botTerm))

if __name__ == '__main__':
    N=7
    x = np.arange(0.01,3.5,0.01)
    y = np.zeros(shape=(2,5,len(x)))
    hoJ = np.arange(0.0,2.5,.5)
    for j in range(len(hoJ)):
        print hoJ[j]
        H_hat = get_H_hat(N,hoJ[j])
        w,v = np.linalg.eig(H_hat)
        for T in range(len(x)):
            y[0,j,T] = E(w,x[T],N)
            y[1,j,T] = C(w,x[T],N)
    for i in range(5):
        plt.plot(x,y[0,i,:],label="h/J="+str(i*0.5))
        plt.legend(loc=2)
        plt.xlabel("T/J")
        plt.ylabel("E/J")
    plt.show()
    for i in range(5):
        plt.plot(x,y[1,i,:],label="h/J="+str(i*0.5))
        plt.legend(loc=1)
        plt.xlabel("T/J")
        plt.ylabel("C/J")
    plt.show()
