from ngsolve import *
from netgen.geom2d import unit_square

ngsglobals.msg_level = 1

mesh = Mesh( unit_square.GenerateMesh())

V = H1(mesh)

m = BilinearForm (V, symmetric=True)
m.Assemble()

mstar = m.mat.CreateMatrix()
mstarinv = mstar.Inverse(V.FreeDofs(),inverse="sparsecholesky")

u = GridFunction (V)
update = u.vec.CreateVector()

u.vec.data = mstarinv * update # RuntimeError: NGSolve exception: FlatArray::operator[]: index 5 out of range [0,4]
