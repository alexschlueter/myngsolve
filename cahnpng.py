from netgen.geom2d import unit_square, SplineGeometry
from ngsolve import *
from ngsapps.utils import *
import random

# convert \( -delay 14 `j=0; for i in $(ls aa.*.png); do if [ $(( $j%2 )) -eq 0 ]; then echo $i\[20%\]; fi; j=$(( $j+1)); done` \) \( -delay 800 'aaend.png[20%]' \) anim2.gif

order = 3

png_width = 400
png_height = 200

maxdim = max(png_width, png_height)
xscale = png_width / maxdim
yscale = png_height / maxdim
print(xscale,yscale)

initial_roughness = 30.0

tau = 1e-6
tend = 3

lamdba = 1e-2
M = 1

vtkoutput = True

def sqr(x):
    return x * x

from PIL import Image

def set_initial_conditions(result_gridfunc, mesh):
    # c0 = GridFunction(result_gridfunc.space)
    c0 = CoefficientFunction(0.0)
    total_mass = 0.0
    # vec_storage = c0.vec.CreateVector()
    # vec_storage[:] = 0.0

    print("setting initial conditions")
    for i in range(80):
        # print("\rtotal mass = {:10.6e}".format(total_mass), end="")
        center_x = xscale * random.random()
        center_y = yscale * random.random()
        thinness_x = initial_roughness * (1+random.random())
        thinness_y = initial_roughness * (1+random.random())
        c0 += exp(-(sqr(thinness_x) * sqr(x-center_x) + sqr(thinness_y) * sqr(y-center_y)))
        # vec_storage.data += c0.vec
        # c0.vec.data = vec_storage

        # cut off above 1.0
    result_gridfunc.Set(IfPos(c0-0.5,0.5,c0))
    # total_mass = Integrate(result_gridfunc,mesh,VOL)

    print()

# add gaussians with random positions and widths until we reach total mass >= 0.5
def set_ic_from_png():
    im = Image.open("ngsolve.png")
    pix = im.load()

    total_mass = 0.0
    func = CoefficientFunction(0.0)

    print("setting initial conditions")

    for i in range(0,png_width):
        for j in range(0,png_height):
            val = pix[i,j][0]/255
            # print (i,"-",j,"-",val,"-",pix[i,j])
            if val == 1.0:
                continue
            center_x = xscale * (i+0.5)/png_width
            center_y = yscale * (1.0 - (j+0.5)/png_height)
            thinness = 130.0
            func += 0.125*exp(-(sqr(thinness) * sqr(x-center_x) + sqr(thinness) * sqr(y-center_y)))

    # result_gridfunc.Set(IfPos(func-1.0,1.0,func))
    return func

# mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
# mesh = Mesh(unit_square.GenerateMesh(maxh=0.03))
geo = SplineGeometry()
geo.AddRectangle((0,0), (xscale,yscale))
mesh = Mesh(geo.GenerateMesh(maxh=0.03))

V = H1(mesh, order=order)
fes = FESpace([V, V])
c, mu = fes.TrialFunction()
q, v = fes.TestFunction()

a = BilinearForm(fes)
a += SymbolicBFI(tau * M * grad(mu) * grad(q))
a += SymbolicBFI(mu * v)
a += SymbolicBFI(-200 * (c - 3 * c * c + 2 * c * c * c) * v)
a += SymbolicBFI(-lamdba * grad(c) * grad(v))

b = BilinearForm(fes)
b += SymbolicBFI(c * q)

b.Assemble()

mstar = b.mat.CreateMatrix()

s = GridFunction(fes)
vec_storage = s.components[0].vec.CreateVector()

# set_initial_conditions(s.components[0],mesh)
# func = RandomCF(0.0,0.2) + set_ic_from_png()
func = set_ic_from_png()
s.components[0].Set(IfPos(func - 1.0, 1.0, func))
Draw(s.components[0], mesh, "c")
# vec_storage.data = s.components[0].vec.data
# set_ic_from_png(s.components[0])
# s.components[0].vec.data += vec_storage.data
# Redraw()
#alternative: use random number generator pointwise:
# s.components[0].Set(GeneratePythonCF(lambda x,y,z: x))
s.components[1].Set(CoefficientFunction(0.0))

rhs = s.vec.CreateVector()
sold = s.vec.CreateVector()
As = s.vec.CreateVector()
w = s.vec.CreateVector()

# Draw(s.components[1], mesh, "mu")
# Draw(s.components[0], mesh, "c")

if vtkoutput:
    vtk = VTKOutput(ma=mesh,coefs=[s.components[0]],names=["c"],filename="cahnhilliard_",subdivision=4)
    vtk.Do()

input("Press any key...")
# implicit Euler
t = 0.0
while t < tend:
    print("\n\nt = {:10.6e}".format(t))

    sold.data = s.vec
    wnorm = 1e99

    # Newton solver
    while wnorm > 1e-9:
        rhs.data = b.mat * sold
        rhs.data -= b.mat * s.vec
        a.Apply(s.vec,As)
        rhs.data -= As
        a.AssembleLinearization(s.vec)

        mstar.AsVector().data = b.mat.AsVector() + a.mat.AsVector()
        invmat = mstar.Inverse()
        w.data = invmat * rhs
        wnorm = w.Norm()
        print("|w| = {:7.3e} ".format(wnorm),end="")
        s.vec.data += w

    t += tau
    Redraw(blocking=False)
    if vtkoutput:
        vtk.Do()
