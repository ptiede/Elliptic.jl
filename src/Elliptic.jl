module Elliptic
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

include("jacobi.jl")
include("slatec.jl")

function E(phi::Number, m::Number)
    if isnan(phi) || isnan(m) return NaN end
    if m < zero(m) || m > one(m) throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*E(m) - _E(cos(mod(phi2,pi)), m)
    end
    _E(sin(phi), m)
end
function _E(sinphi::Number, m::Number)
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
#E(phi::Real, m::Real) = E(Number(phi), Number(m))


"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
function ellipke(m::Number)
  if m < one(m)
        y = 1 - m
        drf,ierr1 = SLATEC.DRF(0, y, 1)
        drd,ierr2 = SLATEC.DRD(0, y, 1)
        @assert ierr1 == 0 && ierr2 == 0
        return drf, drf - m*drd/3
    elseif m == one(m)
      return Inf, one(m)
    elseif isnan(m)
        return NaN, NaN
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end

E(m::Number) = last(ellipke(m))
#E(x::Float32) = Float32(E(Number(x)))
#E(x::Real) = E(Number(x))

# assumes 0 ≤ m ≤ 1
function rawF(sinphi::Number, m::Number)
  if abs(sinphi) == one(sinphi) && m == one(sinphi) return sign(sinphi)*Inf end
    sinphi2 = sinphi^2
    drf,ierr = SLATEC.DRF(1 - sinphi2, 1 - m*sinphi2, 1)
    @assert ierr == zero(ierr)
    sinphi*drf
end

function F(phi::Number, m::Number)
    if isnan(phi) || isnan(m) return NaN end
    if m < zero(m) || m > one(m) throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        # Abramowitz & Stegun (17.4.3)
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*K(m) - rawF(cos(mod(phi2,pi)), m)
    end
    rawF(sin(phi), m)
end
#F(phi::Real, m::Real) = F(Number(phi), Number(m))

function K(m::Number)
    if m < 1.
        drf,ierr = SLATEC.DRF(0, 1 - m, 1)
        @assert ierr == 0
        return drf
    elseif m == one(m)
        return Inf
    elseif isnan(m)
        return NaN
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end
#K(x::Float32) = Float32(K(Number(x)))
#K(x::Real) = K(Number(x))

function Pi(n::Number, phi::Number, m::Number)
    if isnan(n) || isnan(phi) || isnan(m) return NaN end
    if m < zero(m) || m > one(m) throw(DomainError(m, "argument m not in [0,1]")) end
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1 - sinphi2
    y = 1 - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1)
    drj,ierr2 = SLATEC.DRJ(cosphi2, y, 1, 1 - n*sinphi2)
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
#Pi(n::Real, phi::Real, m::Real) = Pi(Number(n), Number(phi), Number(m))
Π = Pi

function ellipj(u::Number, m::Number, tol::Number)
    phi = Jacobi.am(u, m, tol)
    s = sin(phi)
    c = cos(phi)
    d = sqrt(1. - m*s^2)
    s, c, d
end
ellipj(u::Number, m::Number) = ellipj(u, m, eps())
#ellipj(u::Real, m::Real) = ellipj(Number(u), Number(m))

end # module
