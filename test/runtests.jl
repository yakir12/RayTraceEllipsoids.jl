using Test
using RayTraceEllipsoids, LinearAlgebra


using Unitful
import Unitful:m, °, μm
# Light properties
source_distance = 1m
source_angle    = 17°

# Eye morphology
aperture        = 251μm
morphing_factor = 1.12

## Volumes
### Refractive indices

ri.water                         = 1.334
ri.cornea                        = 1.37
ri.lens                          = 1.42
ri.lens2distal_retina            = 1.35
ri.distal_retina2proximal_retina = 1.35
ri.gap                           = 1.34

### Thicknesses

thick.cornea                        = 23μm
thick.lens                          = 215μm
thick.lens2distal_retina            = 6μm
thick.distal_retina2proximal_retina = 81μm
thick.gap                           = 49μm

## Interfaces
### Radii of the ellipsoids that form the interface

radii.cornea_distal   = (244μm, 227μm)
radii.lens_distal     = (195μm, 213μm)
radii.lens_proximal   = 337μm
radii.distal_retina   = 337μm
radii.proximal_retina = 337μm
radii.mirror = 417μm

r = Ray(V3c(1., 2, 10), V3c(0., 0, -1), 1.)
s = Ellipsoid(V3c(1., 2, 2), V3c(5., 4, 2), 1.)
c = OpticUnit(s, true, 1.4, RayTraceEllipsoid.NonRetina(0.1))
RayTraceEllipsoid.raytrace!(r, c)


r = Ray(V3c(1., 2, 10), V3c(0., 0, -1), 1.)
cs = V3c{Float64}(1,2,2)
rs = V3c{Float64}(5,4,2)
s = Ellipsoid(cs, rs, 1.)
c = OpticUnit(s, true, 1.4, RayTraceEllipsoid.Retina(cs, rs, 0.01, 0.1))
function fun()
r = Ray(V3c(1., 2, 10), V3c(0., 0, -1), 1.)
RayTraceEllipsoid.raytrace!(r, c)
end

    l, p = RayTraceEllipsoid.distancepoint(r, c.surface)

    RayTraceEllipsoid.absorb!(r, l, c.medium)

    r.orig = p

    bend!(r, c.interface)

#=using BenchmarkTools, RayTraceEllipsoid, UnitfulAngles
s = Ellipsoid(V3c(1., 2, 2), V3c(5., 4, 2), V3c(0., 0, 1), pi/2)
ou = OpticUnit(s, false, 1., true, "a", ac = 0.05)
function fun()
    r = Ray(V3c(1., 2, 10), V3c(0., 0, -1), 1.)
    for i in 1:7
        a = raytrace!(r, ou)
    end
end
@btime fun() # 271.179 ns (1 allocation: 64 bytes)

s = Ellipsoid(V3c{BigFloat}(1, 2, 2), V3c{BigFloat}(5., 4, 2), V3c{BigFloat}(0., 0, 1), BigFloat(pi)/2)
ou = OpticUnit(s, false, BigFloat(1), true, "a", ac = BigFloat(0.05))
function fun()
    r = Ray(V3c{BigFloat}(1., 2, 10), V3c{BigFloat}(0., 0, -1), BigFloat(1))
    for i in 1:7
        a = raytrace!(r, ou)
    end
end
@btime fun() # 271.179 ns (1 allocation: 64 bytes)
=#

#=using Makie, RayTraceEllipsoid, LinearAlgebra

s = Ellipsoid(V64(1, 2, 2), V64(5, 4, 2), V64(0, 0, 1), pi/2)
n = 100
θs = range(0, stop = 1, length = n)
φs = range(0, stop = 2, length = n)
xs = zeros(n,n)
ys = zeros(n,n)
zs = zeros(n,n)
for i in 1:n, j in 1:n
    θ = θs[i]
    φ = φs[j]
    x = cospi(φ)*sinpi(θ)
    y = sinpi(φ)*sinpi(θ)
    z = cospi(θ)
    pt = V64(x, y, z)
    pt = s.uncenter_unscale(pt)
    xs[i,j], ys[i,j], zs[i,j] = pt
end
surface(xs, ys, zs, color = fill(RGBAf0(1, 0, 0, 0.5),n,n ))#color = :red, transparency = true,  alpha = .5)
oc = V64(4,4,10)
r = Ray(oc, V64(0, 0, -1), 1.)
l = advance!(r, s)
points = [
     Point3f0(oc...) => Point3f0(r.orig)
    ]
linesegments!(points, color = :blue, linewidth = 2)=#






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
