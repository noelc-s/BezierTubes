import sys
sys.path.insert(0, '../build')

import numpy as np
import bezier
import matplotlib.pyplot as plt

b = bezier.Bezier(3,2,1)
A_u = np.zeros((32, 4))
b_u = np.zeros(32)
b.input_constraints(A_u,b_u,[ 1.02006, 0.425478],[0,0],[0])
x0 = np.array([0,0])
 
B = b.hyp2vert(A_u[:,0:2],b_u-A_u[:,2:]@x0)
B = np.vstack((B, B[0,:]))

F = b.hyp2vert(A_u[:,2:4],b_u-A_u[:,0:2]@x0)
F = np.vstack((F, F[0,:]))

t = np.linspace(0,1)

plt.figure()
xs = F[:,0]
ys = F[:,1]
plt.scatter(xs,ys) 
xs = B[:,0]
ys = B[:,1]
plt.scatter(xs,ys) 
for v in B:
  print(np.hstack((v,x0)))
  b_t = b.B(t, np.linalg.inv(b.D.T)@np.hstack((v,x0)))
  plt.plot(b_t[:,0],b_t[:,1],'b')
for v in F:
  b_t = b.B(t, np.linalg.inv(b.D.T)@np.hstack((x0,v)))
  plt.plot(b_t[:,0],b_t[:,1],'g')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show() # if you need...

