module Jacobi

export am,
    sn, cn, dn, nn,
    sd, cd, dd, nd,
    sc, cc, dc, nc,
    ss, cs, ds, ns

# Abramowitz & Stegun, section 16.4, p571
const _ambuf = Array{Float64}(undef, 10)
function _am(u::Number, m::Number, tol::Number)
    if u == 0 return 0. end

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4
        return u - 1/4*m*(u - 1/2*sin(2*u))
    end
    m1 = 1 - m
    if m1 < sqrt_tol
        # A&S 16.15.4
        t = tanh(u)
        return asin(t) + 1/4*m1*(t - u*(1 - t^2))*cosh(u)
    end

    a,b,c,n = 1, sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        @assert n < 10
        a,b,c,n = 1/2*(a+b), sqrt(a*b), 0.5*(a-b), n+1
        _ambuf[n] = c/a
    end

    phi = ldexp(a*u, n)
    for i = n:-1:1
        phi = 0.5*(phi + asin(_ambuf[i]*sin(phi)))
    end
    phi
end
_am(u::T, m::T) where {T} = _am(u, m, eps(T))

"""
    am(u::Real, m::Real, [tol::Real=eps(Float64)])

Returns amplitude, φ, such that u = F(φ | m)

Landen sequence with convergence to `tol` used if `√(tol) ≤ m ≤ 1 - √(tol)`
"""
function am(u::Number, m::Number, tol::Number)
    (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
    _am(u, m, tol)
end
am(u::T, m::T) where {T} = am(u, m, eps(T))

for (f,a,b,c) in ((:sn, :(sin(phi)),                :(sqrtmu1*s), :(sqrt(mu)*sin(phi))),
                  (:cn, :(cos(phi)),                :(cos(phi)),  :(sqrt(1 - mu*sin(phi)^2))),
                  (:dn, :(sqrt(1 - m*sin(phi)^2)), :(1),        :(cos(phi))))
    @eval begin
        function ($f)(u, m)
            # Abramowitz & Stegun, section 16.10, p573
            lt0 = m < 0.
            gt1 = m > 1.
            if !(lt0 || gt1)
                phi = _am(u,m)
                return $a
            elseif lt0
                mu1 = inv(1 - m)
                mu = -m*mu1
                sqrtmu1 = sqrt(mu1)
                v = u/sqrtmu1
                phi = _am(v,mu)
                s = sin(phi)
                return ($b)/sqrt(1 - mu*s^2)
            elseif gt1
                mu = 1/m
                v = u*sqrt(m)
                phi = _am(v,mu)
                return $c
            end
        end
    end
end

xn = ((:s,:(sn(u,m))), (:c,:(cn(u,m))), (:d,:(dn(u,m))), (:n,:(1.)))
for (p,num) in xn, (q,den) in xn
    f = Symbol(p, q)
    @eval begin
        """
            $($f)(u::Real, m::Real)

        Compute the Jacobi elliptic function $($f)(u | m)
        """
        ($f)(u::Real, m::Real) = ($f)(u, m)
    end

    if (p == q)
        @eval ($f)(::T, ::T) where {T} = 1.0
    elseif (q != :n)
        @eval ($f)(u::T, m::T) where {T} = ($num)/($den)
    end
end

end # module
