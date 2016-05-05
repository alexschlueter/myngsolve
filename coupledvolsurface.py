# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.geom2d import MakeCircle
from netgen.geom2d import SplineGeometry
from netgen.csg import unit_cube

ngsglobals.msg_level = 1

dt = 0.001
tend = 1.0

dL = 1.0 # 0.01
dl = 1.0 # 0.02
gamma = 1.0 #2
lamdba = 1.0 #4

# generate a triangular mesh of mesh-size 0.2
# mesh = Mesh (unit_square.GenerateMesh(maxh=0.2))
circ = SplineGeometry()
MakeCircle(circ, (0,0), 1, bc=1)
# mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2))
mesh = Mesh(circ.GenerateMesh(maxh=0.2))
bnd = [i for i, n in enumerate(mesh.GetBoundaries())]
print(str(bnd))

# H1-conforming finite element space
V = H1(mesh, order=3)
W = H1(mesh, order=3, flags={"definedon": [], "definedonbound": [1,2,3,4,5,6], "dirichlet": [1,2,3,4,5,6]})



fes = FESpace([V, W])

L,l = fes.TrialFunction()
v,w = fes.TestFunction()

# the bilinear-form
a = BilinearForm (fes, symmetric=True)
a += SymbolicBFI(dL * L.Deriv() * v.Deriv())
a += SymbolicBFI(dl * l.Trace().Deriv() * w.Trace().Deriv(), BND)
a += SymbolicBFI((lamdba * L - gamma * l.Trace()) * (v - w.Trace()), BND)

c = BilinearForm(fes, symmetric=True)
c += SymbolicBFI(L * v)
c += SymbolicBFI(l.Trace() * w.Trace(), BND)

a.Assemble()
print("a geht")
c.Assemble()
print("c")

# sL = GridFunction(V)
# sL.Set(0.5 * (x * x + y * y))
# sl = GridFunction(W)
# sl.Set(0.5 * (1 + x))
s = GridFunction(fes)
s.components[0].Set(0.5 * (x * x + y * y))
s.components[1].Set(0.5 * (1 + x), boundary=True)

mstar = a.mat.CreateMatrix()
mstar.AsVector().data = c.mat.AsVector() + dt * (a.mat.AsVector())
alldofs = BitArray(fes.ndof)
for i in range(fes.ndof):
    alldofs.Set(i)

invmat = mstar.Inverse(alldofs)
rhs = s.vec.CreateVector()
Draw(s.components[1], mesh, "l")
Draw(s.components[0], mesh, "L")
t = 0.0
input("")
while t < tend:
    rhs.data = c.mat * s.vec
    s.vec.data = invmat * rhs
    t += dt
    Redraw(blocking=True)
    input("")


# # rhs = u.vec.CreateVector()

# # Draw (u)
# # t=0.0
# # while t < tend:
# #     rhs.data = dt * f.vec 
# #     rhs.data += m.mat * u.vec
# #     u.vec.data = invmat * rhs
# #     t += dt
# #     Redraw()
# #     input("")
