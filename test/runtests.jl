using Test
using RayTraceEllipsoids, LinearAlgebra, CoordinateTransformations

@test all(1:10^5) do _
    oz, l, dz = rand(4) .- 0.5
    positive = rand(Bool)
    a = RayTraceEllipsoids.isgoodz(oz, l, dz, positive)
    z = oz + l*dz
    b = z > 0 ? positive : !positive
    a == b
end


@test all(1:10^5) do _
    z = rand()
    orig = V3(0.,0,z + 1)
    dir = V3(0.,0,-1)
    pos = true
    RayTraceEllipsoids.distance(orig, dir, pos) == z
end

@test_throws RayTraceEllipsoids.DeadRay() RayTraceEllipsoids.distance(V3(2,2,2.), V3(0,0,-1.), true)

@test all(1:10^5) do _
    z = rand()
    c = V3(1.,1,1)
    r = V3(1.,2,3)
    h = true
    s = Ellipsoid(c, r, h)
    orig = V3(1,1,c[3] + r[3] + z)
    dir = V3(0.,0,-1)
    r = Ray(orig, dir)
    RayTraceEllipsoids.distance(r, s) ≈ z
end


@test all(1:10^5) do _
    z = rand()
    c = V3(1.0,1,1)
    r = V3(1.0,2,3)
    h = true
    s = Ellipsoid(c, r, h)
    orig = V3(1.0,1,c[3] + r[3] + z)
    dir = V3(0.0,0,-1)
    r = Ray(orig, dir)
    RayTraceEllipsoids.distance(r, s) ≈ z
end



@test all(1:10^5) do _
    orig = V3(0.,0,0)
    l = rand()
    dir = V3(0.,0,-1)
    o = V3(rand(3)...)
    ucs = Translation(o)
    rorig = o
    RayTraceEllipsoids._getl(orig, l, dir, ucs, rorig) == l
end

#1,3
r = Ray(V3(2.,2,2), V3(0.,0,-1))
s = Cylinder(1.5,-2.,2.)
L = 4.
@test RayTraceEllipsoids.incylinder(r, s, L) == 0
#2
r = Ray(V3(1.,1,2), V3(0.,0,-1))
s = Cylinder(1.5,-2.,2.)
L = 4.
@test RayTraceEllipsoids.incylinder(r, s, L) == 4
#4
r = Ray(V3(0.,0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == sqrt(2)
#5
r = Ray(V3(1.,-1,2), V3(0.,1,-1))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == 0
#6
r = Ray(V3(-4tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) ≈ sqrt(2)
#7
r = Ray(V3(-2tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) ≈ sqrt(2^2+2^2)
#8
r = Ray(V3(-6tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
s = Cylinder(1.,-2.,2.)
L = sqrt(4^2+4^2)
@test RayTraceEllipsoids.incylinder(r, s, L) == 0


s = Cylinder(1.,-2.,2.)
m = Retina(s, 0.1)
r = Ray(V3(-2tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
L = sqrt(2^2+2^2)
signal = (1 - RayTraceEllipsoids.absorption(L, 0.1))
RayTraceEllipsoids.register!(r, l, m)
@test m.signal[] == signal


s = Cylinder(1.,-2.,2.)
m = Retina(s, 0.1)
r = Ray(V3(-2tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
o2 = r.orig + r.dir*l
RayTraceEllipsoids.absorbmove!(r, l, m)
i = RayTraceEllipsoids.absorption(l, 0.1)
@test r.int == i
@test r.orig == o2

s = Cylinder(1.,-2.,2.)
m = Retina(s, 1 - RayTraceEllipsoids.ENERGYTHRESHOLD)
r = Ray(V3(-2tand(45),0,2), V3(sqrt(1/2),0,-sqrt(1/2)))
l = sqrt(4^2+4^2)
o2 = r.orig + r.dir*l
@test_throws RayTraceEllipsoids.DeadRay() RayTraceEllipsoids.absorbmove!(r, l, m)

@test all(1:10^5) do _
    z = rand()
    c = V3(1.0,1,1)
    r = V3(2.0,2,2)
    orig = V3(c[1] + 1, c[2] + 1, c[3] + r[3] + z)
    dir = V3(0.0,0,-1)
    p1 = copy(dir)
    ray = Ray(orig, dir)
    h = true
    ss = Ellipsoid(c, r, h)
    l = RayTraceEllipsoids.distance(ray, ss)
    s = Cylinder(r[1], c[2] - r[2], c[3] + r[3])
    m = Retina(s, 0.5)
    RayTraceEllipsoids.absorbmove!(ray, l, m)
    ou = OpticUnit(ss, false, 0.9, m)
    RayTraceEllipsoids.bend!(ray, ou.interface)
    p2 = ray.dir
    acosd(dot(p1, p2) / sqrt(norm(p1)*norm(p2))) ≈ 45 - asind(sind(45)*0.9)
end

@test all(1:10^5) do _
    d = 10*rand()
    r = rand()
    cz = rand()
    d += cz
    l = Light(d, 2r, r, r, cz, 0)
    ray = l(1)
    θ = atan(r/d)
    ray.dir ≈ V3(sin(θ), 0, -cos(θ)) && ray.orig == V3(0, 0, d + cz)
end

@test all(1:10^5) do _
    d = 10*rand()
    r = rand()
    cz = rand()
    d += cz
    l = Light(d, 2r, r, r, cz, 0)
    ray = l(1, 0.5)
    ray.dir ≈ V3(0, 0, -1) && ray.orig == V3(0, 0, d + cz)
end


@test all(1:10^5) do _
    z = rand()
    c = V3(1.0,1,1)
    r = V3(2.0,2,2)
    orig = V3(c[1] + 1, c[2] + 1, c[3] + r[3] + z)
    dir = V3(0.0,0,-1)
    p1 = copy(dir)
    ray = Ray(orig, dir)
    h = true
    ss = Ellipsoid(c, r, h)
    l = RayTraceEllipsoids.distance(ray, ss)
    s = Cylinder(r[1], c[2] - r[2], c[3] + r[3])
    m = Retina(s, 0.5)
    ou = OpticUnit(ss, false, 0.9, m)
    raytrace!(ray, ou)
    p2 = ray.dir
    acosd(dot(p1, p2) / sqrt(norm(p1)*norm(p2))) ≈ 45 - asind(sind(45)*0.9) && norm(ray.orig - orig) == l
end

