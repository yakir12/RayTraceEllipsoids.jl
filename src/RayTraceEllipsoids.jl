module RayTraceEllipsoid

using CoordinateTransformations, StaticArrays, LinearAlgebra

export V3c, Ray, Ellipsoid, Interface, OpticUnit, raytrace!, Medium, NonRetina, Retina, Interface, Glued, Active, Near, Far, Light, Cylinder, cspheroid

const ENERGYTHRESHOLD = 0.01

"""
    V3c = SVector{3,Float64}

Points and directions are just a point in 3D space best described as a static V3ctor, `SVector`. `V3c` is an alias for that.
"""
const V3c = SVector{3}

"""
    Ray(orig::V3c, dir::V3c)

The main Ray type with a ray origin, `orig`, and direction, `dir`. The ray's direction gets normalized.
"""
mutable struct Ray{T <: AbstractFloat}
    orig::V3c{T}
    dir::V3c{T}
    int::T
    Ray(o::V3c{T}, d::V3c{T}) where {T <: AbstractFloat} = new{T}(o, normalize(d), one(T))
end

"""
    Ellipsoid(c::V3c, r::V3c, dir::V3c, open::Float64)

    An ellipsoid with a center, `c`, and radii, `r`, as well as a direction (gets automatically normalized), `dir`, and an opening angle, `open`, creating a dome (or window), that the ellipsoid is defined in. Note that `open` is the angle between the dome's edge and the direction of the dome (so actually half the opening angle) and is defined in **some angular units** (using UnitfulAngles, for example: u"°").

`Ellipsoid` has 6 additional fields all relating to various spatial transformations that convert the ellipsoid to a unit-sphere and back again. These are all `CoordinateTransformations`.

# Examples
For an ellipsoid upwards-pointing hemisphere with a center at (1,2,3), and radii (4,5,6):
```jldoctest
julia> Ellipsoid(V3c(1,2,3), V3c(4,5,6), V3c(0,0,1), 0.0)
Rayden.Ellipsoid([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [0.0, 0.0, 1.0], 0.0, Translation(-1.0, -2.0, -3.0), LinearMap([0.25 0.0 0.0; 0.0 0.2 0.0; 0.0 0.0 0.166667]), AffineMap([0.25 0.0 0.0; 0.0 0.2 0.0; 0.0 0.0 0.166667], [-0.25, -0.4, -0.5]), Translation(1.0, 2.0, 3.0), LinearMap([4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0]), AffineMap([4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0], [1.0, 2.0, 3.0]))
```
"""
struct Ellipsoid{T <: AbstractFloat}
    c::V3c{T}
    r::V3c{T}
    h::T # the point on the vertical axis from which the dome starts
    name::String
    # all of the following are transformations
    hᵀ::T # the point on the vertical axis from which the dome starts
    center::Translation{V3c{T}} # translate the ellipsoid to zero
    scale::LinearMap{SDiagonal{3,T}}
    center_scale::AffineMap{SDiagonal{3,T},V3c{T}} # translate and scale to a unit-sphere
    # uncenter::Translation{V3c{T}}
    # unscale::LinearMap{SDiagonal{3,T}}
    uncenter_unscale::AffineMap{SDiagonal{3,T},V3c{T}}
    function Ellipsoid(c::V3c{T}, r::V3c{T}, h::T, name::String) where {T <: AbstractFloat}
        uncenter = Translation(c)
        unscale = LinearMap(SDiagonal(r))
        center = inv(uncenter)
        scale = inv(unscale)
        center_scale = scale∘center
        uncenter_unscale = inv(center_scale)
        hᵀ = h/r[3]
        new{T}(c, r, h, name, hᵀ, center, scale, center_scale, uncenter_unscale)
    end
end

cspheroid(cz, rxy, rz, h, name) = Ellipsoid(V3c(zero(cz), zero(cz), cz), V3c(rxy, rxy, rz), h, name)

function isgoodz(oz, l, dz, hᵀ)
    z = oz + l*dz
    if hᵀ > 0
        z ≥ hᵀ
    else
        z ≤ hᵀ
    end
end

struct DeadRay <: Exception end

"""
    distance(orig::V3c, dir::V3c)

Return the two distances between a point with origin, `orig`, and direction, `dir`, and the two (potentially identical) intersection points with a unit-sphere. 
"""
function distance(orig::V3c{T}, dir::V3c{T}, hᵀ::T) where {T <: AbstractFloat}
    b = -orig⋅dir
    disc = b^2 - orig⋅orig + 1
    if disc ≥ 0
        d = sqrt(disc)
        t2 = b + d
        if t2 ≥ 0
            t1 = b - d
            if t1 > 0 
                isgoodz(orig[3], t1, dir[3], hᵀ) && return t1
                isgoodz(orig[3], t2, dir[3], hᵀ) && return t2
            else
                isgoodz(orig[3], t2, dir[3], hᵀ) && return t2
            end
        end
    end
    throw(DeadRay())
end

function distance(r::Ray{T}, s::Ellipsoid{T}) where {T <: AbstractFloat}
    orig = s.center_scale(r.orig) # transform the ray's origin according to the ellipsoid
    dir = normalize(s.scale(r.dir)) # scale the ray's direction as well
    l = distance(orig, dir, s.hᵀ)
    o = orig + l*dir
    p = s.uncenter_unscale(o)
    norm(p - r.orig)
end

struct Cylinder{T <: AbstractFloat}
    scale::LinearMap{SDiagonal{3,T}}
    center_scale::AffineMap{SDiagonal{3,T},V3c{T}}
    uncenter_unscale::AffineMap{SDiagonal{3,T},V3c{T}}
    function Cylinder(rxy::T, z1::T, z2::T) where {T <: AbstractFloat}
        cz = (z1 + z2)/2
        uncenter = Translation(V3c(0, 0, cz))
        @assert z1 < z2 "cylinder z1 seems bigger than z2"
        rz = (z2 - z1)/2
        unscale = LinearMap(SDiagonal(rxy, rxy, rz))
        center = inv(uncenter)
        scale = inv(unscale)
        center_scale = scale∘center
        uncenter_unscale = inv(center_scale)
        new{T}(scale, center_scale, uncenter_unscale)
    end
end


function _getl(orig, l, dir, ucs, rorig)
    pᵗ = orig + l*dir
    p2 = ucs(pᵗ)
    norm(p2 - rorig)
end

function incylinder(r::Ray{T}, s::Cylinder{T}, L::T) where {T <: AbstractFloat}
    orig = s.center_scale(r.orig)
    dir = normalize(s.scale(r.dir))
    a = dir[1]^2 + dir[2]^2
    if a == 0
        if orig[1]^2 + orig[2]^2 < 1
            return L
        else
            return zero(T)
        end
    end
    b = 2orig[1]*dir[1] + 2orig[2]*dir[2]
    c = orig[1]^2 + orig[2]^2 - 1
    Δ = b^2 - 4a*c
    if Δ ≤ 0
        return zero(T)
    end
    sqrtΔ = sqrt(Δ)
    l1 = (-b - sqrtΔ)/2a
    l2 = (-b + sqrtΔ)/2a
    l1, l2 = l1 < l2 ? (l1, l2) : (l2, l1)
    if l2 < 0
        return zero(T)
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
            return zero(T)
        end
        l2 = _getl(orig, l2, dir, s.uncenter_unscale, r.orig)
        if l2 > L
            return L - l1
        else
            return l2 - l1
        end
    end
end

abstract type Medium{T <: AbstractFloat} end

struct NonRetina{T <: AbstractFloat} <: Medium{T}
    absorption_coefficient::T
    name::String
end

struct Retina{T <: AbstractFloat} <: Medium{T}
    cylinder::Cylinder{T}
    absorption_coefficient::T
    signal::Base.RefValue{T}
    n::Base.RefValue{Int}
    name::String
    function Retina(cylinder::Cylinder{T}, absorption_coefficient::T, name::String) where {T <: AbstractFloat}
        new{T}(cylinder, absorption_coefficient, Ref(zero(T)), Ref(0), name)
    end
end

register!(r, l, m::NonRetina) = nothing

absorption(l, absorption_coefficient) = exp(-l*absorption_coefficient)

function register!(r, l, m::Retina)
    l = incylinder(r, m.cylinder, l)
    if !iszero(l)
        m.signal[] += r.int*(1 - absorption(l, m.absorption_coefficient))
        m.n[] += 1
    end
    nothing
end

function absorbmove!(r, l, m)
    register!(r, l, m)
    r.int *= absorption(l, m.absorption_coefficient)
    r.orig += l*r.dir
end



#=struct Retina{T <: AbstractFloat} <: Medium{T}
    uncenter_unscale::AffineMap{SDiagonal{3,T},V3c{T}}
    center_scale::AffineMap{SDiagonal{3,T},V3c{T}} # translate and scale to a unit-sphere
    σ::WeightedVariance{T}
    step::T
    absorbed::T
    CCD::OffsetArray{T,2,Array{T,2}}
    n::Base.RefValue{Int64}
    function Retina(uncenter_unscale::AffineMap{SDiagonal{3,T},V3c{T}}, step::T, absorption_coefficient::T) where {T <: AbstractFloat}
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

abstract type Interface{T <: AbstractFloat} end

struct Glued{T <: AbstractFloat} <: Interface{T}
end

"""
    Interface(normal::AffineMap, n::Float64)

Build an optical interface from a AffineMap that transforms a point on the ellipsoid of the interface to the normal at that point, `normal`, and the refractive index ratio between the inside and the outside of the ellipsoid, `n`.
"""
struct Active{T <: AbstractFloat} <: Interface{T}
    normal::AffineMap{SDiagonal{3, T}, V3c{T}}
    n::T
    n²::T
    Active(normal::AffineMap{SDiagonal{3,T},V3c{T}}, n::T) where {T <: AbstractFloat} = new{T}(normal, n, n^2)
end

refract(n, dir, a, b, N) = n*dir + (n*a - sqrt(1 - b))*N

reflect(dir, a, N) = dir + 2a*N

"""
    bend!(r::Ray, i::Interface)

Refract or reflect a ray with an interface. Update the direction of the ray. Returns if event failed, which is always `true` (see `raytrace!` for details why that is so).
"""
function bend!(r::Ray{T}, i::Active{T}) where {T <: AbstractFloat}
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

bend!(r::Ray{T}, i::Glued{T}) where {T <: AbstractFloat} = nothing

"""
    OpticUnit(surface::Ellipsoid, pointin::Bool, n::Float64, register::Bool, name::String)

Build an optical unit from a surface. `pointin` dictates if the normal should be pointing in or out. `n` is the refractive index. `register` indicates whether this unit should register intersection points (thus functions as a retina).  
"""
struct OpticUnit{T <: AbstractFloat}
    surface::Ellipsoid{T}
    medium::Medium{T}
    interface::Interface{T}

    function OpticUnit(surface::Ellipsoid{T}, pointin::Bool, n::T, medium::Medium{T}) where {T <: AbstractFloat}
        interface = if n ≠ 1
            i = (-one(T))^pointin
            dir = LinearMap(SDiagonal(i, i, i))
            normal = dir∘surface.scale∘surface.scale∘surface.center
            Active(normal, n)
        else
            Glued{T}()
        end
        new{T}(surface, medium, interface)
    end
end

"""
    raytrace!(r::Ray, c::OpticUnit)

Advance a ray to the intersection point with an ellipsoid. If the intersection was successful, bend the ray. Updates the ray accordingly. Returns intersection failure.
"""
function raytrace!(r::Ray{T}, c::OpticUnit{T}) where {T <: AbstractFloat}
    l = distance(r, c.surface)
    absorbmove!(r, l, c.medium)
    r.int < ENERGYTHRESHOLD && throw(DeadRay())
    bend!(r, c.interface)
end

function raytrace!(r::Ray{T}, ous::Vector{OpticUnit{T}}) where {T <: AbstractFloat}
    try
        foreach(ous) do ou
            raytrace!(r, ou)
        end
    catch ex
        ex isa DeadRay || throw(ex)
    end
end

# light

abstract type Light{T} end

struct Near{T <: AbstractFloat} <: Light{T}
    cosa::T
    orig::V3c{T}
    rotm::LinearMap{RotY{T}}
end

struct Far{T <: AbstractFloat} <: Light{T}
    R::T
    z::T
    dir::V3c{T}
    rotm::AffineMap{RotY{T}, V3c{T}}
end

function Light(distance::T, aperture::T, rxy::T, rz::T, cz::T, θ::T) where {T <: AbstractFloat}
    r = aperture/2
    h = sqrt((1 - r^2/rxy^2)*rz^2) + cz
    rotm = LinearMap(RotY(θ))
    tran = recenter(rotm, V3c{T}(0,0,h))
    if distance < 1e6 # todo, check accuracy here
        orig = tran(V3c{T}(0,0,distance))
        p2 = V3c{T}(r,0,h) - orig 
        cosalpha1 = dot(normalize(-orig), normalize(p2))
        p2 = V3c{T}(-r,0,h) - orig
        cosalpha2 = dot(normalize(-orig), normalize(p2))
        cosalpha3 = cos(atan(r/(distance - h)))
        cosa = minimum([cosalpha1, cosalpha2, cosalpha3])
        Near{T}(cosa, orig, rotm)
    else
        Far{T}(r, cz + rz + 100, rotm(V3c{T}(0,0,-1)), tran)
    end
end

function (l::Near{T})(b1::T, b2::T) where {T <: AbstractFloat}
    z = b1*(l.cosa - 1) - l.cosa
    k = sqrt(1 - z^2)
    theta = 2b2
    dir = l.rotm(V3c{T}(k*cospi(theta), k*sinpi(theta), z))
    Ray(l.orig, dir)
end

(l::Near{T})() where {T <: AbstractFloat} = l(rand(T), rand(T))

function (l::Far{T})(b1::T, b2::T) where {T <: AbstractFloat}
    k = l.R*sqrt(b1)
    theta = 2b2
    orig = l.rotm(V3c{T}(k*cospi(theta), k*sinpi(theta), l.z))
    Ray(orig, l.dir)
end

(l::Far{T})() where {T <: AbstractFloat} = l(rand(T), rand(T))

#=struct ScallopEye{T <: AbstractFloat, L <: Light}
    aperture::T
    l::L
    ous::Vector{OpticUnit{T}}
    n::Int

    function OpticalSystem(aperture::T, distance::T, θ::T, n::Int, tissues::Vector{Tissue}, membranes::Vector{Membrane}, morphz::T) where {T <: AbstractFloat}
        new{T, L}(aperture, l, ous, n)
    end
end=#

function raytrace!(l::Light{T}, ous::Vector{OpticUnit{T}}, n::Int) where {T <: AbstractFloat}
    foreach(1:n) do _
        r = l()
        raytrace!(r, ous)
    end
end


end #module
