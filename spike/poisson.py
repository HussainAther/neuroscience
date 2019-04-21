from numpy import *
from numpy.matlib import *
from numpy.linalg import *
from matplotlib.pyplot import *
from matplotlib.delaunay.triangulate import Triangulation
from math import factorial

"""
Poisson process between vectors A and b to see if there is some event that occurs in both.

Poisson's equation is obtained from adding a source term to the right-hand-side of Laplace's equation:

∂2p/∂x2+∂2p/∂y2 = b

So, unlike the Laplace equation, there is some finite value inside the field that affects the solution.
Poisson's equation acts to "relax" the initial sources in the field.
"""

def meshgrid(xs, ys, npoints):
    """
    Meshgrids are 2-D grid coordinates for the simulation.
    """
    # randomly choose some points
    rng = random.RandomState(1234567890)
    rx = rng.uniform(xs[0], xs[1], size=npoints)
    ry = rng.uniform(ys[0], ys[1], size=npoints)
    # only take points in domain
    nx, ny = [], []
    for x,y in zip(rx,ry):
        if in_domain(x,y):
            nx.append(x)
            ny.append(y)
    # Delaunay triangulation
    tri = Triangulation(array(nx), array(ny))
    return tri

def A_e(v):
    # take vertices of element and return contribution to A
    G = vstack((ones((1,3)), v.T)).I * vstack((zeros((1,2)),eye(2)))
    return det(vstack((ones((1,3)), v.T))) * G * G.T / 2

def b_e(v):
    # take vertices of element and return contribution to b
    vS = v.sum(axis=0)/3.0 # Centre of gravity
    return f(vS) * ((v[1,0]-v[0,0])*(v[2,1]-v[0,1])-(v[2,0]-v[0,0])*(v[1,1]-v[0,1])) / 6.0

def poisson(tri, boundary):
    # get elements and vertices from meshgrid
    elements = tri.triangle_nodes
    vertices = vstack((tri.x,tri.y)).T
    # number of vertices and elements
    N = vertices.shape[0]
    E = elements.shape[0]
    #Loop over elements and assemble LHS and RHS
    A = zeros((N,N))
    b = zeros((N,1))
    for j in range(E):
        index = (elements[j,:]).tolist()
        A[ix_(index,index)] += A_e(vertices[index,:])
        b[index] += b_e(vertices[index,:])
    # find the 'free' vertices that we need to solve for
    free = list(set(range(len(vertices))) - set(boundary))
    # initialise solution to zero so 'non-free' vertices are by default zero
    u = zeros((N,1))
    # solve for 'free' vertices.
    u[free] = solve(A[ix_(free,free)],b[free])
    return array(u)

def f(v):
    # the RHS f
    x, y = v
    f = 2.0*cos(10.0*x)*sin(10.0*y) + sin(10.0*x*y)
    return 1

def in_domain(x,y):
    # is a point in the domain?
    return sqrt(x**2 + y**2) <= 1

xs = (-1.,1.)
ys = (-1.,1.)
npoints = 1000

# generate meshgrid and determine boundary vertices
tri = meshgrid(xs, ys, npoints)
boundary = tri.hull

# solve Poisson equation
u  = poisson(tri, boundary).flatten()

# interpolate values and plot a nice image
lpi = tri.linear_interpolator(u)
z = lpi[ys[0]:ys[1]:complex(0,npoints),
        xs[0]:xs[1]:complex(0,npoints)]
z = where(isinf(z), 0.0, z)
extent = (xs[0], xs[1], ys[0], ys[1])
ioff()
clf()
imshow(nan_to_num(z), interpolation='bilinear', extent=extent, origin='lower')
#show()
savefig("sol.png", bb_inches="tight")

"""
Compare the Poisson model to actual data using Fano factors, interspike interval distributions,
and coefficients of variation.

Fano factor describes relationship between mean spike count over a given interval and the spike-
count variance.
"""

a = [ 1.00000000000000000000, 0.57721566490153286061, -0.65587807152025388108,
         -0.04200263503409523553, 0.16653861138229148950, -0.04219773455554433675,
         -0.00962197152787697356, 0.00721894324666309954, -0.00116516759185906511,
         -0.00021524167411495097, 0.00012805028238811619, -0.00002013485478078824,
         -0.00000125049348214267, 0.00000113302723198170, -0.00000020563384169776,
          0.00000000611609510448, 0.00000000500200764447, -0.00000000118127457049,
          0.00000000010434267117, 0.00000000000778226344, -0.00000000000369680562,
          0.00000000000051003703, -0.00000000000002058326, -0.00000000000000534812,
          0.00000000000000122678, -0.00000000000000011813, 0.00000000000000000119,
          0.00000000000000000141, -0.00000000000000000023, 0.00000000000000000002]

def gamma(n, product=1):
    """
    Gamma function for some number n. Based off distribution.
    """
    y = float(n) - 1.0;
    s = a[-1];
    for an in a[-2::-1]:
        s = s * y + an
    return 1.0 / sm

def interspikeInterval(x, alpha, beta):
    """
    Gamma distribution is better for modeling interspike intrevals than the exponential Poisson
    distribution. The refractoriness makes short spike intervals less liekly than the
    Poisson model would predict.
    
    We use a x is gamma-distributed with a shape alpha and rate beta.
    In this case, beta is the time constant, x is the rate of firing, and alpha
    is the refractory exponential.
    """
    num = (beta**alpha) * (x**(alpha-1)) * np.exp(-beta*x) # numerator
    den = gamma(alpha) # denominator
    return num/den
