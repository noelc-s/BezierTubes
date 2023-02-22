import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import sys
sys.path.append('/home/amber/repos/shortest-paths-in-graphs-of-convex-sets')
from spp.convex_sets import Singleton, Polyhedron, Ellipsoid
from spp.convex_functions import TwoNorm, SquaredTwoNorm
from spp.graph import GraphOfConvexSets
from spp.shortest_path import ShortestPathProblem
from pydrake.all import MosekSolver
import os
os.environ["MOSEKLM_LICENSE_FILE"] = "/home/amber/mosek/mosek.lic"

exec(open("poly.txt").read())

MosekSolver.AcquireLicense()

# convex sets
singletons = (
  Singleton((-1.9, -1.9)),
  Singleton((1.9, 1.9)),
)
sets = singletons + polyhedra

# label for the vertices
vertices = ['s', 't']
vertices += [f'p{i}' for i in range(len(polyhedra))]

# add convex sets to the graph
G = GraphOfConvexSets()
G.add_sets(sets, vertices)
G.set_source('s')
G.set_target('t')

# edges
H = np.hstack((np.eye(2), -np.eye(2)))
l = TwoNorm(H)
for u, vs in edges.items():
  print(u)
  print(": ")
  for v in vs:
    print(v)
    G.add_edge(u, v, l)
  print("----")

spp = ShortestPathProblem(G, relaxation=0)
sol = spp.solve()



path = np.concatenate((sol.primal.x[0], sol.primal.x[0]),axis=0)
print(path)
EdgeList = np.array([[0, 0]])
for k, phi in enumerate(sol.primal.phi):
            if phi > 1 - 1e-3:
                edge = [G.vertices.index(vertex) for vertex in G.edges[k]]
                print(sol.primal.x[edge])
                print(np.linalg.norm(sol.primal.x[edge][0] - sol.primal.x[edge][1]))
                if (np.linalg.norm(sol.primal.x[edge][0] - sol.primal.x[edge][1]) > 1e-3):
                  EdgeList = np.vstack([EdgeList, edge])
                  path = np.vstack([path, sol.primal.x[edge].reshape(4)])
path = path[1:]
EdgeList = EdgeList[1:]
print(path)
print(EdgeList)

with open('gcs_optimized_path.m', 'w') as f:
  f.write("path="+np.array2string(path))
  f.write(';\n')
  f.write("edgeTraversal="+np.array2string(EdgeList)+";")
