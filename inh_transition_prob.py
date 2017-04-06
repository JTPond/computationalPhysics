import numpy as np
from matplotlib import pyplot as plt

def prob(x,y):
    if abs(y-x) == 1 and y != 0 and y != 11:
        return .5*((2.*x - 10. -1.)/(10.0))**2.0
    elif x==y and (y != 1 and y != 10):
        return 1.0 -((2.*x - 10. -1.)/(10.0))**2.0
    elif x==y and (y == 1 or y == 10):
        return 1.0 -.5*((2.*x - 10. -1.)/(10.0))**2.0 
    else:
        return 0.0

P = np.ndarray(shape=(10,10))

for i in range(10):
    for j in range(10):
        P[i,j] = prob(j+1,i+1)

print "\n",P,"\n"

w,v = np.linalg.eig(P)

print w,"\n" #[ 0.20131988  0.20131988  0.68932593  0.68932773  1.          0.99219895  0.97267342  0.97472871  0.89462368  0.89448183]

#w is a vector of eigan values, index 4 = 1.
#The numpy.linalg documentation states the w[i] is the e-value
#   of e-state v[:,i] don't ask me why... but it uses LAPACK to do the calculation

print v[:,4],"\n" #[ 0.00866704  0.01432715  0.02808121  0.07800337  0.7020303   0.7020303  0.07800337  0.02808121  0.01432715  0.00866704]

print (v[:,4]**2),"\n" #[  7.51175942e-05   2.05267195e-04   7.88554457e-04   6.08452513e-03  4.92846536e-01   4.92846536e-01   6.08452513e-03   7.88554457e-04  2.05267195e-04   7.51175942e-05]
print (v[:,4]**2).sum(),"\n" #1.0

plt.plot(range(1,11,1),v[:,4]**2,linestyle='dashed',marker='*') #Plots the distrebution... 
plt.show()
