using Test
using RayTraceEllipsoids, LinearAlgebra, CoordinateTransformations

@test all(1:10^5) do _
    oz, l, dz, hᵀ = rand(4) .- 0.5
    a = RayTraceEllipsoids.isgoodz(oz, l, dz, hᵀ) 
    z = oz + l*dz
    b = hᵀ > 0 ? z ≥ hᵀ : z ≤ hᵀ
    a == b
end


@test all(1:10^5) do _
    z = rand()
    orig = V3c(0.,0,z + 1)
    dir = V3c(0.,0,-1)
    hᵀ = 0.5
    RayTraceEllipsoids.distance(orig, dir, hᵀ) == z
end

@test_throws RayTraceEllipsoids.DeadRay() RayTraceEllipsoids.distance(V3c(2,2,2.), V3c(0,0,-1.), 0.5)

@test all(1:10^5) do _
    z = rand()
    c = V3c(1.,1,1)
    r = V3c(1.,2,3)
    h = .5
    name = "a"
    s = Ellipsoid(c, r, h, name)
    orig = V3c(1,1,c[3] + r[3] + z)
    dir = V3c(0.,0,-1)
    r = Ray(orig, dir)
    RayTraceEllipsoids.distance(r, s) ≈ z
end

@test all(1:10^5) do _
    orig = V3c(0.,0,0)
    l = rand()
    dir = V3c(0.,0,-1)
    o = V3c(rand(3)...)
    ucs = Translation(o)
    rorig = o
    RayTraceEllipsoids._getl(orig, l, dir, ucs, rorig) == l
end

#1,3
r = Ray(V3c(2.,2,2), V3c(0.,0,-1))
s = Cylinder(1.5,-2.,2.)
L = 4.
@test RayTraceEllipsoids.incylinder(r, s, L) == 0
#2
r = Ray(V3c(1.,1,2), V3c(0.,0,-1))
s = Cylinder(1.5,-2.,2.)
L = 4.
@test RayTraceEllipsoids.incylinder(r, s, L) == 4
#4
r = Ray(V3c(0.,0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == sqrt(2)
#5
r = Ray(V3c(1.,-1,2), V3c(0.,1,-1))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == 0
#6
r = Ray(V3c(-4tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) ≈ sqrt(2)
#7
r = Ray(V3c(-2tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) ≈ sqrt(2^2+2^2)
#8
r = Ray(V3c(-6tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == 0


s = Cylinder(1.,-2.,2.)
m = Retina(s, 0.1, "a")
r = Ray(V3c(-2tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
L = sqrt(2^2+2^2)
signal = (1 - RayTraceEllipsoids.absorption(L, 0.1))
RayTraceEllipsoids.register!(r, l, m)
@test m.signal[] == signal
@test m.n[] == 1


s = Cylinder(1.,-2.,2.)
m = Retina(s, 0.1, "a")
r = Ray(V3c(-2tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
o2 = r.orig + r.dir*l
RayTraceEllipsoids.absorbmove!(r, l, m)
i = RayTraceEllipsoids.absorption(l, 0.1)
@test r.int == i
@test r.orig == o2

s = Cylinder(1.,-2.,2.)
m = Retina(s, 1 - RayTraceEllipsoids.ENERGYTHRESHOLD, "a")
r = Ray(V3c(-2tand(45),0,2), V3c(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
o2 = r.orig + r.dir*l
@test_throws RayTraceEllipsoids.DeadRay() RayTraceEllipsoids.absorbmove!(r, l, m)

