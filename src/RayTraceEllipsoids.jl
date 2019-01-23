module RayTraceEllipsoids

using CoordinateTransformations
using StaticArrays: SVector, SDiagonal
using LinearAlgebra: normalize, ⋅, norm
# import Optim
# using OffsetArrays
# using WeightedOnlineStats
# using AngleBetweenVectors

export V3, Ray, EmptyRay, TraceRay, DeadRay, AbstractRay, Ellipsoid, Interface, OpticUnit, raytrace!, trace!, Medium, NonRetina, Retina, Interface, Glued, Active, Light, Cylinder, cspheroid, Signal

const ENERGYTHRESHOLD = 0.01

"""
    V3 = SVector{3,Float64}

Points and directions are just a point in 3D space best described as a static `SVector`. `V3` is an alias for that.
"""
const V3 = SVector{3, Float64}

abstract type AbstractRay end

struct DeadRay <: AbstractRay end

"""
    Ray(orig::V3, dir::V3, int::Float64)

The main Ray type with a ray origin, `orig`, a direction, `dir`, and an intensity `int`. The ray's direction gets normalized.
"""
struct Ray <: AbstractRay
    orig::V3
    dir::V3
    int::Float64
end
Ray(o, d) = Ray(o, d, 1.0)
# Ray() = Ray(zero(V3), zero(V3))

"""
    EmptyRay(orig::V3, dir::V3, int::Float64)

A Ray type for plotting. Same as `Ray` but doesn't get registered in retinas. Mainly used for plotting.
"""
struct EmptyRay <: AbstractRay
    orig::V3
    dir::V3
    int::Float64
end
EmptyRay(o, d) = EmptyRay(o, d, 1.0)


"""
    TraceRay(orig::V3, dir::V3, int::Float64)

A Ray type for plotting. Same as `Ray` but doesn't get registered in retinas. Mainly used for plotting.
"""
struct TraceRay <: AbstractRay
    orig::V3
    dir::V3
end

# TraceRay() = TraceRay(zero(V3), zero(V3))

# struct PSFRay <: AbstractRay
    # orig::V3
    # dir::V3
    # int::Float64
    # PSFRay(o::V3, d::V3, int::Float64) = new(o, d, int)
# end
# PSFRay() = PSFRay(zero(V3), zero(V3))

"""
    Ellipsoid(c::V3, r::V3, dir::V3, open::Float64)

    An ellipsoid with a center, `c`, and radii, `r`, as well as a direction (gets automatically normalized), `dir`, and an opening angle, `open`, creating a dome (or window), that the ellipsoid is defined in. Note that `open` is the angle between the dome's edge and the direction of the dome (so actually half the opening angle) and is defined in **some angular units** (using UnitfulAngles, for example: u"°").

`Ellipsoid` has 6 additional fields all relating to various spatial transformations that convert the ellipsoid to a unit-sphere and back again. These are all `CoordinateTransformations`.
```
"""

struct Ellipsoid
    positive::Bool # is the (normalized) dome in the positive end of the z-axis?
    center::Translation{V3} # translate the ellipsoid to zero
    scale::LinearMap{SDiagonal{3,Float64}} # scale to a unit-sphere
    center_scale::AffineMap{SDiagonal{3,Float64},V3} # translate and scale to a unit-sphere
    uncenter_unscale::AffineMap{SDiagonal{3,Float64},V3}

    function Ellipsoid(c::V3, r::V3, positive::Bool)
        uncenter = Translation(c)
        unscale = LinearMap(SDiagonal(r))
        center = inv(uncenter)
        scale = inv(unscale)
        center_scale = scale∘center
        uncenter_unscale = inv(center_scale)
        new(positive, center, scale, center_scale, uncenter_unscale)
    end

end

cspheroid(cz::Real, rxy::Real, rz::Real, h::Bool) = Ellipsoid(V3(0.0, 0.0, cz), V3(rxy, rxy, rz), h)

isgoodz(oz, l, dz, positive::Bool) = oz > -l*dz ? positive : !positive

"""
    distance(orig::V3, dir::V3)

Return the two distances between a point with origin, `orig`, and direction, `dir`, and the two (potentially identical) intersection points with a unit-sphere. 
"""
function distance(orig::V3, dir::V3, positive::Bool)
    b = -orig⋅dir
    disc = b^2 - orig⋅orig + 1
    if disc ≥ 0
        d = sqrt(disc)
        t2 = b + d
        if t2 ≥ 0
            t1 = b - d
            if t1 > 0 
                isgoodz(orig[3], t1, dir[3], positive) && return t1
                isgoodz(orig[3], t2, dir[3], positive) && return t2
            else                                                   
                isgoodz(orig[3], t2, dir[3], positive) && return t2
            end
        end
    end
    nothing
end

function distance(r::AbstractRay, s::Ellipsoid)
    orig = s.center_scale(r.orig) # transform the ray's origin according to the ellipsoid
    dir = normalize(s.scale(r.dir)) # scale the ray's direction as well
    l = distance(orig, dir, s.positive)
    untransform(l, orig, dir, r.orig, s.uncenter_unscale)
end
untransform(l::Nothing, orig, dir, rorig, uncenter_unscale) = l
function untransform(l::Float64, orig, dir, rorig, uncenter_unscale)
    o = orig + l*dir
    p = uncenter_unscale(o)
    norm(p - rorig)
end

struct Cylinder
    scale::LinearMap{SDiagonal{3,Float64}}
    center_scale::AffineMap{SDiagonal{3,Float64},V3}
    uncenter_unscale::AffineMap{SDiagonal{3,Float64},V3}

    function Cylinder(rxy::Real, z1::Real, z2::Real)
        cz = (z1 + z2)/2
        uncenter = Translation(V3(0.0, 0.0, cz))
        @assert z1 < z2 "cylinder z1 seems bigger than z2"
        rz = (z2 - z1)/2
        unscale = LinearMap(SDiagonal(rxy, rxy, rz))
        center = inv(uncenter)
        scale = inv(unscale)
        center_scale = scale∘center
        uncenter_unscale = inv(center_scale)
        new(scale, center_scale, uncenter_unscale)
    end

end


function _getl(orig::V3, l::Float64, dir::V3, ucs, rorig::V3)
    pᵗ = orig + l*dir
    p2 = ucs(pᵗ)
    norm(p2 - rorig)
end

function incylinder(r::Ray, s::Cylinder, L::Float64)
    orig = s.center_scale(r.orig)
    dir = normalize(s.scale(r.dir))
    a = dir[1]^2 + dir[2]^2
    if a == 0
        if orig[1]^2 + orig[2]^2 < 1
            return L
        else
            return 0.0
        end
    end
    b = 2orig[1]*dir[1] + 2orig[2]*dir[2]
    c = orig[1]^2 + orig[2]^2 - 1
    Δ = b^2 - 4a*c
    if Δ ≤ 0
        return 0.0
    end
    sqrtΔ = sqrt(Δ)
    l1 = (-b - sqrtΔ)/2a
    l2 = (-b + sqrtΔ)/2a
    l1, l2 = l1 < l2 ? (l1, l2) : (l2, l1)
    if l2 < 0
        return 0.0
    end
    if l1 < 0
        l = _getl(orig, l2, dir, s.uncenter_unscale, r.orig)
        if 0 < l ≤ L
            return l
        else
            return L
        end
    else
        l1 = _getl(orig, l1, dir, s.uncenter_unscale, r.orig)
        if l1 > L
            return 0.0
        end
        l2 = _getl(orig, l2, dir, s.uncenter_unscale, r.orig)
        if l2 > L
            return L - l1
        else
            return l2 - l1
        end
    end
end

abstract type Medium end

struct NonRetina <: Medium
    absorption_coefficient::Float64
end

mutable struct Signal
    photoreceptor::Float64
    retina::Float64
    Signal() = new(0.0, 0.0)
end

struct Retina <: Medium
    cylinder::Cylinder
    absorption_coefficient::Float64
    signal::Signal
    # photons::Base.RefValue{Int}
    # pixels::OffsetArray{Float64}
    # n::Base.RefValue{Int}
    # var::WeightedMean{Float64}
    # σ²::WeightedVariance{Float64}
    # data::Vector{Float64}
    # weight::Vector{Float64}
    # nodalz::Float64
    function Retina(cylinder::Cylinder, absorption_coefficient::Float64)
        # new(absorption_coefficient, WeightedMean{Float64}(), Float64[], Float64[], nodalz)
        # pixels = OffsetArray{Float64}(undef, -400:400, -400:400)
        # fill!(pixels, 0.0)
        new(cylinder, absorption_coefficient, Signal())#, WeightedVariance{Float64}())
    end
end

const LOG10EXP = 1/(log10(exp(1)))

function absorption(l::Float64, absorption_coefficient::Float64) 
    kl = absorption_coefficient*l
    kl/(LOG10EXP + kl)
end
# absorption(l::Float64, absorption_coefficient::Float64) = exp(-l*absorption_coefficient)

function register!(m::Retina, r::Ray, l, w)
    linside = incylinder(r, m.cylinder, l)
    m.signal.photoreceptor += r.int*absorption(linside, m.absorption_coefficient)
    m.signal.retina += w
    # if linside > 0
    #     m.photons[] += 1
    # end
    # fit!(m.σ², r.orig[1], w/2)
    # fit!(m.σ², r.orig[2], w/2)
    nothing
end

register!(m, r, l, w) = nothing

abstract type Interface end

struct Glued <: Interface
    passthrough::Float64
end

"""
    Interface(normal::AffineMap, n::Float64)

Build an optical interface from a AffineMap that transforms a point on the ellipsoid of the interface to the normal at that point, `normal`, and the refractive index ratio between the inside and the outside of the ellipsoid, `n`.
"""
struct Active <: Interface
    passthrough::Float64
    normal::AffineMap{SDiagonal{3, Float64}, V3}
    n::Float64
    n²::Float64
    Active(pass::Float64, normal::AffineMap{SDiagonal{3, Float64}, V3}, n::Float64) = new(pass, normal, n, n^2)
end

refract(n, dir, a, b, N) = n*dir + (n*a - sqrt(1 - b))*N

reflect(dir, a, N) = dir + 2a*N

"""
    newdirection(r::Ray, i::Interface)

Refract or reflect a ray with an interface. Update the direction of the ray. Returns if event failed, which is always `true` (see `raytrace!` for details why that is so).
"""
function newdirection(orig, dir, i::Active)
    N = normalize(i.normal(orig))
    a = -dir⋅N
    b = i.n²*(1 - a^2)
    dir = if b ≤ 1
        refract(i.n, dir, a, b, N)
    else
        reflect(dir, a, N)
    end
    normalize(dir)
end

newdirection(orig, dir, i::Glued) = dir

"""
    OpticUnit(surface::Ellipsoid, pointin::Bool, n::Float64, register::Bool, name::String)

Build an optical unit from a surface. `pointin` dictates if the normal should be pointing in or out. `n` is the refractive index. `register` indicates whether this unit should register intersection points (thus functions as a retina).  
"""
struct OpticUnit{T <: Medium}
    surface::Ellipsoid
    medium::T
    interface::Interface

    function OpticUnit(surface::Ellipsoid, pointin::Bool, n::Float64, medium::T, passthrough::Float64) where {T<:Medium}
        interface = if n ≠ 1
            i = Float64((-1)^pointin)
            dir = LinearMap(SDiagonal(i, i, i))
            normal = dir∘surface.scale∘surface.scale∘surface.center
            Active(passthrough, normal, n)
        else
            Glued(passthrough)
        end
        new{T}(surface, medium, interface)
    end
end

trace!(c::OpticUnit, r::DeadRay) = r

"""
    trace!(c::OpticUnit, r::Ray)

Advance a ray to the intersection point with an ellipsoid. If the intersection was successful, bend the ray. Updates the ray accordingly. Returns intersection failure.
"""
function trace!(c::OpticUnit, r::AbstractRay)
    l = distance(r, c.surface)
    mutate(l, r, c)
end

newlocation(orig, l, dir) = orig + l*dir

mutate(l::Nothing, r, c) = DeadRay()
function mutate(l::Float64, r::TraceRay, c)
    orig = newlocation(r.orig, l, r.dir)
    dir = newdirection(orig, r.dir, c.interface)
    TraceRay(orig, dir)
end
function mutate(l::Float64, r::AbstractRay, c)
    w = absorbed(r, l, c.medium)
    int = newintensity(r, w)
    _mutate(int, w, r, l, c)
end
_mutate(int::Nothing, w, r, l, c) = DeadRay()
function _mutate(int, w, r, l, c)
    register!(c.medium, r, l, w)
    orig = newlocation(r.orig, l, r.dir)
    int = passed(c.interface.passthrough, int)
    __mutate(int, orig, r, c.interface)
end

function __mutate(int::Float64, orig, r, interface)
    dir = newdirection(orig, r.dir, interface)
    typeof(r)(orig, dir, int)
end

__mutate(int::Nothing, orig, r, interface) = DeadRay()

function passed(passthrough, int)
    int *= passthrough
    int < ENERGYTHRESHOLD ?  nothing : int
end

absorbed(r, l, m) = r.int*absorption(l, m.absorption_coefficient)

function newintensity(r, w) 
    int = r.int - w
    int < ENERGYTHRESHOLD ?  nothing : int
end

raytrace!(ous, r::DeadRay) = r
function raytrace!(ous, r::AbstractRay)
    @inbounds for ou in ous
        r = trace!(ou, r)
    end
    r
end

# light

struct Light
    cosa::Float64
    orig::V3
    rotm::LinearMap{RotY{Float64}}
    nom::Float64
    radius²::Float64

end
function Light(distance::Real, aperture::Real, rc, θ::Real, nodalz::Real)
    nodal = V3(0, 0, nodalz)
    radius = aperture/2
    h = rc.rz*sqrt(1 - (radius/rc.rxy)^2) + rc.cz
    vs = (V3(0, radius, h), V3(radius, 0, h), V3(-radius, 0, h))
    rotm = LinearMap(RotY{Float64}(θ))
    tran = recenter(rotm, nodal)
    orig = tran(V3(0, 0, nodalz + distance))
    v1 = normalize(nodal - orig)
    cosa = Inf
    for vi in vs
        _cosa = getcosa(vi, orig, v1)
        if _cosa < cosa
            cosa = _cosa
        end
    end
    nom = h - orig[3]
    Light(cosa, orig, rotm, nom, radius^2)
end
#=function Light(distance::Real, aperture::Real, rc, θ::Real, γ::Real, nodalz::Real)

    =##=nodal = V3(0, 0, nodalz)
    rotm = LinearMap(RotZY{Float64}(γ, θ))
    tran = recenter(rotm, nodal)
    orig = tran(V3(0, 0, nodalz + distance))

    v1 = normalize(nodal - orig)
    radius = aperture/2
    h = rc.rz*sqrt(1 - (radius/rc.rxy)^2) + rc.cz
    res = Optim.optimize(ξ -> ξ2v(ξ, radius, h, orig, v1), 0, π)
    cosa = Optim.minimum(res)
    acosd(cosa)=##=

    nodal = V3(0, 0, nodalz)
    radius = aperture/2
    h = rc.rz*sqrt(1 - (radius/rc.rxy)^2) + rc.cz
    vs = (V3(0, radius, h), V3(radius, 0, h), V3(-radius, 0, h))
    roty = LinearMap(RotY{Float64}(θ))
    trany = recenter(roty, nodal)
    origy = trany(V3(0, 0, nodalz + distance))
    v1 = normalize(nodal - origy)
    # cosa = minimum(getcosa(vi, origy, v1) for vi in vs)
    cosa = Inf
    for vi in vs
        _cosa = getcosa(vi, origy, v1)
        if _cosa < cosa
            cosa = _cosa
        end
    end
    rotm = LinearMap(RotZY{Float64}(γ, θ))
    tran = recenter(rotm, nodal)
    orig = tran(V3(0, 0, nodalz + distance))
    nom = h - orig[3]
    new(cosa, orig, rotm, nom, radius^2)
end=#

function getcosa(vi, origy, v1)
    vin = normalize(vi - origy)
    vin⋅v1
end

#=function ξ2v(ξ, radius, h, orig, v1)
sinξ, cosξ = sincos(ξ)
v = V3(radius*cosξ, radius*sinξ, h) - orig
normalize(v)⋅v1
end=#



#=function Light(distance::Real, aperture::Real, rxy::Real, rz::Real, cz::Real, θ::Real)
radius = aperture/2
h = rz*sqrt(1 - (radius/rxy)^2)
nodal = V3(0, 0, cz + h)
rotm = LinearMap(RotY{Float64}(θ))
tran = recenter(rotm, nodal)
orig = tran(nodal + V3(0, 0,distance))
v1 = nodal - orig
v2 = V3(radius, 0, h + cz) - orig
cosalpha1 = normalize(v1)⋅normalize(v2)
v2 = V3(-radius, 0, h + cz) - orig
cosalpha2 = normalize(v1)⋅normalize(v2)
cosalpha3 = cos(atan(radius/distance))
cosa = min(cosalpha1, cosalpha2, cosalpha3)
z = cz + h
nom = z - orig[3]
Light(cosa, orig, rotm, nom, nodal, radius + 1e-7) # adding some fuzz factor to the radius
end=#


function (l::Light)(::Type{T}, b1::Real, b2::Real) where {T <: AbstractRay}
    z = b1*(l.cosa - 1) - l.cosa
    k = sqrt(1 - z^2)
    theta = 2b2
    dir = l.rotm(V3(k*cospi(theta), k*sinpi(theta), z))
    t = l.nom/dir[3]
    p0 = l.orig + t*dir
    if p0[1]^2 + p0[2]^2 ≤ l.radius²
        T(l.orig, dir)
    else
        DeadRay()
    end
end

function (l::Light)(::Type{T}, b::Real) where {T <: AbstractRay}
    α = acos(l.cosa)
    a = 2α*(b - 0.5)
    dir = l.rotm(normalize(V3(sin(a), 0.0, -cos(a))))
    t = l.nom/dir[3]
    p0 = l.orig + t*dir
    if p0[1]^2 + p0[2]^2 ≤ l.radius²
        T(l.orig, dir)
    else
        DeadRay()
    end
end

(l::Light)(::Type{T}) where {T <: AbstractRay} = l(T, rand(), rand())

function raytrace!(ous, l::Light; n=10^5) # after some testing I found that 10^5 rays is enough for stable stats
    ni = Ref(0)
    @inbounds for _ in 1:n
        r = l(Ray)
        addone!(ni, r)
        raytrace!(ous, r)
    end
    ni[]
end

addone!(ni, r::DeadRay) = nothing
function addone!(ni, r::AbstractRay)
    ni[] += 1
    nothing
end



end #module

