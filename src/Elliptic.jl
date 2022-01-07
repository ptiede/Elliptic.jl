module Elliptic

# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

include("jacobi.jl")
using .Jacobi

include("slatec.jl")

function E(phi, m)
    T = promote_type(typeof(phi), typeof(m))
    if isnan(phi) || isnan(m) return NaN end
    if m < zero(T)  || m > one(T) throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*E(m) - _E(cos(mod(phi2,pi)), m)
    end
    _E(sin(phi), m)
end
function _E(sinphi, m)
    sinphi2 = sinphi^2
    cosphi2 = 1 - sinphi2
    y = 1 - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1)
    drd,ierr2 = SLATEC.DRD(cosphi2, y, 1)
    if ierr1 == ierr2 == 0
        return sinphi*(drf - m*sinphi2*drd/3)
    elseif ierr1 == ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return sinphi
    end
    NaN
end
#E(phi::Real, m::Real) = E(Float64(phi), Float64(m))


"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
function ellipke(m)
    if m < 1
        y = 1 - m
        drf,ierr1 = SLATEC.DRF(0., y, 1.)
        drd,ierr2 = SLATEC.DRD(0., y, 1.)
        @assert ierr1 == 0 && ierr2 == 0
        return drf, drf - m*drd/3
    elseif m == 1.
        return Inf, 1.
    elseif isnan(m)
        return NaN, NaN
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end
#ellipke(x::Real) = ellipke(Float64(x))

E(m) = ellipke(m)[2]
E(x::Float32) = Float32(E(Float64(x)))
#E(x::Real) = E(Float64(x))

# assumes 0 ≤ m ≤ 1
function rawF(sinphi, m)
    if abs(sinphi) == 1. && m == 1. return sign(sinphi)*Inf end
    sinphi2 = sinphi^2
    drf,ierr = SLATEC.DRF(1 - sinphi2, 1 - m*sinphi2, 1.)
    @assert ierr == 0
    sinphi*drf
end

function F(phi, m)
    if m > 1
        ## Abramowitz & Stegum (17.4.15)
        m12 = sqrt(m)
        theta = asin(m12*sin(phi))
        return 1/m12*_F(theta, 1/m)
        #return NaN
    elseif m < 0
        # Abramowitz & Stegum (17.4.17)
        n = -m
        m12 = 1/sqrt(1+n)
        m1m = n/(1+n)
        return m12*K(m1m) - m12*_F(π/2-phi, m1m)
    end
    return _F(phi, m)
end

function _F(phi, m)
    if isnan(phi) || isnan(m) return NaN end
    #if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        # Abramowitz & Stegun (17.4.3)
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*K(m) - rawF(cos(mod(phi2,pi)), m)
    end
    rawF(sin(phi), m)
end
#F(phi::Real, m::Real) = F(Float64(phi), Float64(m))

function K(m)
    if m < 1.
        drf,ierr = SLATEC.DRF(0., 1. - m, 1.)
        @assert ierr == 0
        return drf
    elseif m == 1.
        return Inf
    elseif isnan(m)
        return NaN
    else
        return NaN
        #throw(DomainError(m, "argument m not <= 1"))
    end
end
K(x::Float32) = Float32(K(Float64(x)))
#K(x::Real) = K(Float64(x))

function Pi(n, phi, m)
    if isnan(n) || isnan(phi) || isnan(m) return NaN end
    if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1. - sinphi2
    y = 1. - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1.)
    drj,ierr2 = SLATEC.DRJ(cosphi2, y, 1., 1. - n*sinphi2)
    if ierr1 == 0 && ierr2 == 0
        return sinphi*(drf + n*sinphi2*drj/3)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return Inf
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinphi2 < tol
        return Inf
    end
    NaN
end
#Pi(n::Real, phi::Real, m::Real) = Pi(Float64(n), Float64(phi), Float64(m))
Π = Pi

# function ellipj(u::Float64, m::Float64, tol::Float64)
    # phi = Jacobi.am(u, m, tol)
    # s = sin(phi)
    # c = cos(phi)
    # d = sqrt(1. - m*s^2)
    # s, c, d
# end
# ellipj(u::Float64, m::Float64) = ellipj(u, m, eps(Float64))
# ellipj(u::Real, m::Real) = ellipj(Float64(u), Float64(m))

include("chainrules.jl")

end # module
