from ngsolve import *

from netgen.csg import unit_cube
from netgen.csg import Sphere
from netgen.csg import Pnt
from netgen.csg import CSGeometry
from netgen.geom2d import SplineGeometry
from netgen.geom2d import unit_square

ngsglobals.msg_level = 1

tau = 1e-3
tend = 1

order = 3

geometry = "square"

if geometry == "cube":
    mesh = Mesh( unit_cube.GenerateMesh(maxh=0.1))
elif geometry == "square":
    mesh = Mesh( unit_square.GenerateMesh(maxh=0.1))
elif geometry == "circle":
    geo = SplineGeometry()
    geo.AddCircle ( (0.0, 0.0), r=1, leftdomain=1, rightdomain=0)
    mesh = Mesh( geo.GenerateMesh(maxh=0.2))
    mesh.Curve(order)
elif geometry == "sphere":
    geom = CSGeometry()
    geom.Add (Sphere (Pnt(0,0,0), 1))
    mesh = Mesh( geom.GenerateMesh(maxh=0.2))
    mesh.Curve(order)


# H1-conforming finite element space
Vbulk = H1(mesh, order=order, dirichlet=[])
Vsurface = H1(mesh, order=order, dirichlet=[1,2,3,4,5,6], flags={ "definedon" : [], "definedonbound" : [1,2,3,4,5,6] })
#, flags ={ "definedon" : [], "definedonbound" : [1,2,3,4]})
print ("Vbulk.ndof ", Vbulk.ndof)
print ("Vsurface.ndof ", Vsurface.ndof)

V = FESpace([Vbulk,Vsurface])

u, ub = V.TrialFunction()
v, vb = V.TestFunction()

gradu = u.Deriv()
gradv = v.Deriv()

gradub = ub.Trace().Deriv()
gradvb = vb.Trace().Deriv()

a = BilinearForm (V, symmetric=True)
a += SymbolicBFI (gradu*gradv)
a += SymbolicBFI (gradub*gradvb,BND)
a += SymbolicBFI ((u.Trace()-ub)*(v.Trace()-vb),BND)
a.Assemble(heapsize=10000000)

m = BilinearForm (V, symmetric=True)
m += SymbolicBFI (u*v)
m += SymbolicBFI (ub*vb,BND)
print("noch")
m.Assemble(heapsize=10000000)

mstar = m.mat.CreateMatrix()
mstar.AsVector().data = tau * a.mat.AsVector() + m.mat.AsVector()

# print(a.mat)

# fb = sin(x)
# fv = cos(x)

f = LinearForm (V)
# f += SymbolicLFI (fb*vb,BND)
# f += SymbolicLFI (fv*v)
print("jo")
f.Assemble()
print("a")


# the solution field 
u = GridFunction (V)

# u.components[0].Set(CoefficientFunction(1.0))
# u.components[1].Set(CoefficientFunction(10.0),boundary=True)
u.components[0].Set(0.5*(x*x+y*y))
u.components[1].Set(0.5*(1+x),boundary=True)
print("b")

#u.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="sparsecholesky") * f.vec
alldofs = BitArray(V.ndof)
for i in range(V.ndof):
    alldofs.Set(i)
mstarinv = mstar.Inverse(alldofs,inverse="sparsecholesky")
print("c")

# print (u.vec)

# # plot the solution (netgen-gui only)
Draw (u.components[1])
Draw (u.components[0])
print("d")

rhs = u.vec.CreateVector()
update = u.vec.CreateVector()
print("e")


bulkconc = []
surfconc = []
sumconc = []
times = []
t = 0 
bulkconc.append(Integrate(u.components[0],mesh,VOL))
surfconc.append(Integrate(u.components[1],mesh,BND))
sumconc.append(bulkconc[-1]+surfconc[-1])
times.append(t)
print("f")
while t < tend:
    t += tau
    # rhs.data = -tau * a.mat * u.vec
    # update.data = mstarinv * rhs
    # u.vec.data += update

    print("i")
    update.data = m.mat * u.vec
    print("j")
    u.vec.data = mstarinv * update
    print("k")
    Redraw(blocking=True)
    print("g")
    print ("\r t = {:10.6e}".format(t),end="")
    print("h")
#    Integrate
    bulkconc.append(Integrate(u.components[0],mesh,VOL))
    surfconc.append(Integrate(u.components[1],mesh,BND))
    sumconc.append(bulkconc[-1]+surfconc[-1])
    times.append(t)
    # input("")

import matplotlib.pyplot as plt

plt.figure(0)
# plt.yscale('log')
plt.xlabel("time")
plt.plot(times,bulkconc,'8-')
plt.plot(times,surfconc,'8-')
plt.plot(times,sumconc,'8-')
plt.legend(["total bulk concentration","total surface concentration","total concentration"])
plt.ion()
plt.show()

