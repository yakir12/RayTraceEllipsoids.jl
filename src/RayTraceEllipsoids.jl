module RayTraceEllipsoids

using CoordinateTransformations
using StaticArrays: SVector, SDiagonal
using LinearAlgebra: normalize, ⋅, norm

export V3, Ray, Ellipsoid, Interface, OpticUnit, raytrace!, Medium, NonRetina, Retina, Interface, Glued, Active, Light, Cylinder, cspheroid

const ENERGYTHRESHOLD = 0.01

"""
    V3 = SVector{3,Float64}

Points and directions are just a point in 3D space best described as a static `SVector`. `V3` is an alias for that.
"""
const V3 = SVector{3, Float64}

abstract type AbstractRay end

"""
    Ray(orig::V3, dir::V3, int::Float64)

The main Ray type with a ray origin, `orig`, a direction, `dir`, and an intensity `int`. The ray's direction gets normalized.
"""
mutable struct Ray <: AbstractRay
    orig::V3
    dir::V3
    int::Float64
    Ray(o::V3, d::V3) = new(o, normalize(d), 1.0)
end

"""
    EmptyRay(orig::V3, dir::V3, int::Float64)

A Ray type for plotting. Same as `Ray` but doesn't get registered in retinas. Mainly used for plotting.
"""
mutable struct EmptyRay <: AbstractRay
    orig::V3
    dir::V3
    int::Float64
    EmptyRay(o::V3, d::V3) = new(o, normalize(d), 1.0)
end

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

isgoodz(oz, l, dz, positive::Bool) = oz > -l*dz ? positive : ~positive

struct DeadRay <: Exception end

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
    throw(DeadRay())
end

function distance(r::AbstractRay, s::Ellipsoid)
    orig = s.center_scale(r.orig) # transform the ray's origin according to the ellipsoid
    dir = normalize(s.scale(r.dir)) # scale the ray's direction as well
    l = distance(orig, dir, s.positive)
    o = orig + l*dir
    p = s.uncenter_unscale(o)
    norm(p - r.orig)
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

struct Retina <: Medium
    cylinder::Cylinder
    absorption_coefficient::Float64
    signal::Base.RefValue{Float64}
    # n::Base.RefValue{Int}
    function Retina(cylinder::Cylinder, absorption_coefficient::Float64)
        new(cylinder, absorption_coefficient, Ref(0.0))#, Ref(0))
    end
end

register!(r::AbstractRay, l::Float64, m::NonRetina) = nothing

absorption(l::Float64, absorption_coefficient::Float64) = exp(-l*absorption_coefficient)

function register!(r::Ray, l::Float64, m::Retina)
    l = incylinder(r, m.cylinder, l)
    if !iszero(l)
        m.signal[] += r.int*(1 - absorption(l, m.absorption_coefficient))
        # m.n[] += 1
    end
    nothing
end

register!(r::EmptyRay, l::Float64, m::Retina) = nothing

function absorbmove!(r::AbstractRay, l::Float64, m::Medium)
    register!(r, l, m)
    r.int *= absorption(l, m.absorption_coefficient)
    r.int < ENERGYTHRESHOLD && throw(DeadRay())
    r.orig += l*r.dir
end



#=struct Retina{T <: AbstractFloat} <: Medium{T}
    uncenter_unscale::AffineMap{SDiagonal{3,T},V3{T}}
    center_scale::AffineMap{SDiagonal{3,T},V3{T}} # translate and scale to a unit-sphere
    σ::WeightedVariance{T}
    step::T
    absorbed::T
    CCD::OffsetArray{T,2,Array{T,2}}
    n::Base.RefValue{Int64}
    function Retina(uncenter_unscale::AffineMap{SDiagonal{3,T},V3{T}}, step::T, absorption_coefficient::T) where {T <: AbstractFloat}
        center_scale = inv(uncenter_unscale)
        σ = WeightedVariance(T)
        absorbed = 1 - absorption(step, absorption_coefficient)
        CCD = OffsetArray{T}(undef, -200:200, -200:200)
        CCD .= zero(T)
        new{T}(uncenter_unscale, center_scale, σ, step, absorbed, CCD, Ref(0))
    end
end

function _absorb!(r, V, m)
    r.orig += V
    pᵀ = m.center_scale(r.orig)
    pᴺ = normalize(pᵀ)
    pˢ = m.uncenter_unscale(pᴺ)
    sig = r.int*m.absorbed
    # WeightedOnlineStats.fit!(σ, pˢ[1:2], fill(1/2, 2))
    WeightedOnlineStats.fit!(m.σ, pˢ[1:2], fill(sig/2, 2))
    r.int -= sig
    i = round(Int, pˢ[1])
    j = round(Int, pˢ[2])
    m.CCD[i, j] += sig
end

function absorb!(r, l, m::Retina)
    ls = range(zero(l), stop = l, step = m.step)
    V = m.step*r.dir
    for i in 1:length(ls)
        _absorb!(r, V, m)
    end
    li = l - last(ls)
    if !iszero(li)
        V = li*r.dir
        _absorb!(r, V, m)
    end
    m.n[] += 1
end=#

abstract type Interface end

struct Glued <: Interface
end

"""
    Interface(normal::AffineMap, n::Float64)

Build an optical interface from a AffineMap that transforms a point on the ellipsoid of the interface to the normal at that point, `normal`, and the refractive index ratio between the inside and the outside of the ellipsoid, `n`.
"""
struct Active <: Interface
    normal::AffineMap{SDiagonal{3, Float64}, V3}
    n::Float64
    n²::Float64
    Active(normal::AffineMap{SDiagonal{3, Float64}, V3}, n::Float64) = new(normal, n, n^2)
end

refract(n, dir, a, b, N) = n*dir + (n*a - sqrt(1 - b))*N

reflect(dir, a, N) = dir + 2a*N

"""
    bend!(r::Ray, i::Interface)

Refract or reflect a ray with an interface. Update the direction of the ray. Returns if event failed, which is always `true` (see `raytrace!` for details why that is so).
"""
function bend!(r::AbstractRay, i::Active)
    N = normalize(i.normal(r.orig))
    a = -r.dir⋅N
    b = i.n²*(1 - a^2)
    dir = if b ≤ 1
        refract(i.n, r.dir, a, b, N)
    else
        reflect(r.dir, a, N)
    end
    r.dir = normalize(dir)
    nothing
end

bend!(r::AbstractRay, i::Glued) = nothing

"""
    OpticUnit(surface::Ellipsoid, pointin::Bool, n::Float64, register::Bool, name::String)

Build an optical unit from a surface. `pointin` dictates if the normal should be pointing in or out. `n` is the refractive index. `register` indicates whether this unit should register intersection points (thus functions as a retina).  
"""
struct OpticUnit
    surface::Ellipsoid
    medium::Medium
    interface::Interface

    function OpticUnit(surface::Ellipsoid, pointin::Bool, n::Float64, medium::Medium)
        interface = if n ≠ 1
            i = Float64((-1)^pointin)
            dir = LinearMap(SDiagonal(i, i, i))
            normal = dir∘surface.scale∘surface.scale∘surface.center
            Active(normal, n)
        else
            Glued()
        end
        new(surface, medium, interface)
    end
end

"""
    raytrace!(r::Ray, c::OpticUnit)

Advance a ray to the intersection point with an ellipsoid. If the intersection was successful, bend the ray. Updates the ray accordingly. Returns intersection failure.
"""
function raytrace!(r::AbstractRay, c::OpticUnit)
    l = distance(r, c.surface)
    absorbmove!(r, l, c.medium)
    bend!(r, c.interface)
end

function raytrace!(r::AbstractRay, ous::Vector{OpticUnit})
    try
        foreach(ous) do ou
            raytrace!(r, ou)
        end
    catch ex
        ex isa DeadRay || throw(ex)
    end
end

# light

struct Light
    cosa::Float64
    orig::V3
    rotm::LinearMap{RotY{Float64}}
    nom::Float64
    nodal::V3
    radius::Float64

    function Light(distance::Real, aperture::Real, rxy::Real, rz::Real, cz::Real, θ::Real)
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
        new(cosa, orig, rotm, nom, nodal, radius + 1e-7) # adding some fuzz factor to the radius
    end

end


function (l::Light)(b1, b2; empty = false)
    z = b1*(l.cosa - 1) - l.cosa
    k = sqrt(1 - z^2)
    theta = 2b2
    dir = l.rotm(V3(k*cospi(theta), k*sinpi(theta), z))
    t = l.nom/dir[3]
    p0 = l.orig + t*dir
    norm(p0 - l.nodal) > l.radius && throw(DeadRay())
    empty ? EmptyRay(l.orig, dir) : Ray(l.orig, dir)
end

(l::Light)(; empty = false) = l(rand(), rand(), empty = empty)

function (l::Light)(b; empty = true)
    α = acos(l.cosa)
    a = 2α*(b - 0.5)
    dir = l.rotm(normalize(V3(sin(a), 0.0, -cos(a))))
    t = l.nom/dir[3]
    p0 = l.orig + t*dir
    norm(p0 - l.nodal) > l.radius && throw(DeadRay())
    empty ? EmptyRay(l.orig, dir) : Ray(l.orig, dir)
end

function raytrace!(b, l::Light, ous::Vector{OpticUnit})
    try
        r = l(b)
        raytrace!(r, ous)
    catch ex
        ex isa DeadRay || throw(ex)
    end
end

function raytrace!(b1, b2, l::Light, ous::Vector{OpticUnit})
    try
        r = l(b1, b2)
        raytrace!(r, ous)
    catch ex
        ex isa DeadRay || throw(ex)
    end
end


end #module
