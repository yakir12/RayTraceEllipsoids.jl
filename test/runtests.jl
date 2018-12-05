using Test
using RayTraceEllipsoids, LinearAlgebra

RayTraceEllipsoids.isgoodz(0, 1, 0.5, 0.5+1)

r = Ray()
@test r.dir == V3c(1,0,0)

r = Ray(V64(1, 2, 10), V64(0, 0, -1), 1.)
s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, 1), pi/2)
l = advance!(r, s)
@test norm(s.center(r.orig)) ≈ s.r[3]

r = Ray(V64(1, 2, 10), V64(0, 0, 1), 1.)
l = advance!(r, s)
@test isinf(l)

r = Ray(V64(1, 2, 10), V64(0, 0, -1), 1.)
s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, -1), pi/2)
l = advance!(r, s)
@test s.center(r.orig)[3] ≈ -s.r[3]

r = Ray(V64(1, 2, 10), V64(-.3, .4, -1), 1.)
dir_org = deepcopy(r.dir)
s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, -1), pi/2)
l = advance!(r, s)
ou = OpticUnit(s, true, 1., true, "a", 2e-8)
bend!(r, ou.interface)
@test dir_org ≈ r.dir

r = Ray(V64(1, 2, 10), V64(-.3, .4, -1), 1.)
dir_org = deepcopy(r.dir)
s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, -1), pi/2)
l = advance!(r, s)
ou = OpticUnit(s, true, 1.5, true, "a", 2e-8)
bend!(r, ou.interface)
@test dir_org ≠ r.dir

r = Ray(V64(1, 2, 10), V64(0, 0, -1), 1.)
dir_org = deepcopy(r.dir)
s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, 1), pi/2)
ou = OpticUnit(s, false, 1., true, "a", 0.05)
a = raytrace!(r, ou)
@test norm(s.center(r.orig)) ≈ s.r[3]
@test dir_org ≈ r.dir
