#! /usr/bin/python
import numpy as np
import itertools
from matplotlib import pyplot as plt

def get_states(N):
    out = []
    for i in itertools.combinations_with_replacement([0,1],N):
        for j in itertools.permutations(i,N):
            o = list(j)
            o.reverse()
            if not o in out:
                out.append(o)
    return out #returns list of all possible states of N spins

def get_E(state):
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

def get_off_diag_P(state1,state2,N,T):
    flips=0
    for i in range(len(state1)):
        if not state1[i] == state2[i]:
            flips +=1
    if not flips == 1:
        return 0.0
    else:
        return min([1.0,np.e**(-(get_E(state1) - get_E(state2))/T)])/float(N)

def get_pi_hat(states,N,T):
    ret = np.zeros(shape=(2**N,2**N))
    for i in range(2**N):
        for j in range(2**N):
            if not i == j:
                ret[i,j] = get_off_diag_P(states[i],states[j],N,T)
    for k in range(2**N):
        if ret[:,k].sum() != 1.0:
            ret[k,k] = 1.0 - ret[:,k].sum()
    return ret

def analytical_Sol(T,N):
    Z  =((2.0*np.cosh(1.0/T))**N + (2.0*np.sinh(1.0/T))**N)
    return (np.e**(float(N)/T))/Z

if __name__ == '__main__':
    for N in [10]:
        states = get_states(N)
        T=np.arange(0.15,3.5,0.05)
        Y = []
        y = []
        for t in T:
            pi_hat = get_pi_hat(states,N,t)
            w,v = np.linalg.eig(pi_hat)
            w = list(w)
            o = v[:,w.index(max(w))]
            Y.append(o[-1]/o.sum())
        y = [analytical_Sol(t,N) for t in T]
        plt.plot(T,y,label="Analytical Solution")
        plt.plot(T,Y,linestyle='',marker='*',label="Master Eq")
    plt.title("U_0(Polarized)")
    plt.xlabel("T/J")
    plt.legend()
    plt.show()
    for N in [3,5,10]:
        states = get_states(N)
        T=np.arange(0.15,3.5,0.05)
        Y = []
        for t in T:
            pi_hat = get_pi_hat(states,N,t)
            w,v = np.linalg.eig(pi_hat)
            w = list(w)
            w.sort()
            o = v[:,w.index(max(w))]
            Y.append(-1.0/(np.log(w[-2])))
        plt.plot(T,Y,marker='*',label="N="+str(N))
    plt.title("Tau")
    plt.xlabel("T/J")
    plt.ylim(0.0,300)
    plt.legend()
    plt.show()
