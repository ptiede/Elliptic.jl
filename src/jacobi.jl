module Jacobi

export ellipj

@inline descstep(m) = m/(1+sqrt(1-m))^2

@generated function shrinkm(m,::Val{N}) where {N}
    # [1, Sec 16.12] / https://dlmf.nist.gov/22.7.i
    quote
        f = one(m)
        Base.Cartesian.@nexprs $N i->begin
            k_i = descstep(m)
            m = k_i^2
            f *= 1+k_i
        end
        return (Base.Cartesian.@ntuple $N k), f, m
    end
end

function ellipj_smallm(u,m)
    # [1, Sec 16.13] / https://dlmf.nist.gov/22.10.ii
    sinu, cosu = sincos(u)
    sn = sinu - m*(u-sinu*cosu)*cosu/4
    cn = cosu + m*(u-sinu*cosu)*sinu/4
    dn = 1 - m*sinu^2/2;
    return sn,cn,dn
end

function ellipj_largem(u,m1)
    # [1, Sec 16.15] / https://dlmf.nist.gov/22.10.ii
    sinhu = sinh(u)
    coshu = cosh(u)
    tanhu = sinhu/coshu
    sechu = 1/coshu
    sn = tanhu + m1*(sinhu*coshu-u)*sechu^2/4
    cn = sechu - m1*(sinhu*coshu-u)*tanhu*sechu/4
    dn = sechu + m1*(sinhu*coshu+u)*tanhu*sechu/4
    return sn,cn,dn
end

function ellipj_growm(sn,cn,dn, k)
    # [1, Sec 16.12] / https://dlmf.nist.gov/22.7.i
    for kk in reverse(k)
        sn,cn,dn = (1+kk)*sn/(1+kk*sn^2),
                       cn*dn/(1+kk*sn^2),
                 (1-kk*sn^2)/(1+kk*sn^2)
                 # ^ Use [1, 16.9.1]. Idea taken from [2]
    end
    return sn,cn,dn
end
function ellipj_shrinkm(sn,cn,dn, k::NTuple{N,<:Any}) where {N}
    # [1, Sec 16.14] / https://dlmf.nist.gov/22.7.ii
    for kk in reverse(k)
        sn,cn,dn = (1+kk)*sn*cn/dn,
                 (cn^2-kk*sn^2)/dn, # Use [1, 16.9.1]
                 (cn^2+kk*sn^2)/dn  # Use [1, 16.9.1]
    end
    return sn,cn,dn
end

function ellipj_viasmallm(u,m,::Val{N}) where {N}
    k,f,m = shrinkm(m,Val{N}())
    sn,cn,dn = ellipj_smallm(u/f,m)
    return ellipj_growm(sn,cn,dn,k)
end
function ellipj_vialargem(u,m,::Val{N}) where {N}
    k,f,m1 = shrinkm(1-m,Val{N}())
    sn,cn,dn = ellipj_largem(u/f,m1)
    return ellipj_shrinkm(sn,cn,dn,k)
end


#----------------
# Pick algorithm

function ellipj_dispatch(u,m, ::Val{N}) where {N}
    if abs(m) <= 1 && real(m) <= 0.5
        return ellipj_viasmallm(u,m, Val{N}())
    elseif abs(1-m) <= 1
        return ellipj_vialargem(u,m, Val{N}())
    elseif imag(m) == 0 && real(m) < 0
        # [1, Sec 16.10]
        sn,cn,dn = ellipj_dispatch(u*sqrt(1-m),-m/(1-m), Val{N}())
        return sn/(dn*sqrt(1-m)), cn/dn, 1/dn
    else
        # [1, Sec 16.11]
        sn,cn,dn = ellipj_dispatch(u*sqrt(m),1/m, Val{N}())
        return sn/sqrt(m), dn, cn
    end
end

Base.@pure function nsteps(M)
    mstart = _convfac(M)
    ε = _vareps(M)
    mp = mstart
    i = 0
    while abs(mp) > ε
        mp = descstep(mp)^2
        i += 1
    end
    return i
end

@inline _convfac(::Type{<:Real}) = 0.5
@inline _convfac(::Type{<:Complex}) = 0.5 + sqrt(3)/2im

@inline _vareps(::Type{<:Complex{T}}) where {T} = sqrt(eps(T))
@inline _vareps(T::Type{<:Real}) = sqrt(eps(T))

#Base.@pure nsteps(ε,::Type{<:Real}) = nsteps(0.5,ε) # Guarantees convergence in [-1,0.5]
#Base.@pure nsteps(ε,::Type{<:Complex}) = nsteps(0.5+sqrt(3)/2im,ε) # This is heuristic.
function ellipj_nsteps(u,m)
    # Compute the number of Landen steps required to reach machine precision.
    # For all FloatXX types, this can be done at compile time, while for
    # BigFloat this has to be done at runtime.
    T = promote_type(typeof(u),typeof(m))
    N = nsteps(T)
    return ellipj_dispatch(u,m,Val{N}())::NTuple{3,T}
end


#-----------------------------------
# Type promotion and special values

function ellipj_check(u,m)
    if isfinite(u) && isfinite(m)
        return ellipj_nsteps(u,m)
    else
        T = promote_type(typeof(u),typeof(m))
        return (T(NaN),T(NaN),T(NaN))
    end
end

ellipj(u::Real,m::Real) = ellipj_check(promote(float(u),float(m))...)


function ellipj(u::Complex,m::Real)
    T = promote_type(real.(typeof.(float.((u,m))))...)
    return ellipj_check(convert(Complex{T},u), convert(T,m))
end
ellipj(u,m::Complex) = ellipj_check(promote(float(u),float(m))...)

"""
    ellipj(u,m) -> sn,cn,dn
Jacobi elliptic functions `sn`, `cn` and `dn`.
Convenience function `pq(u,m)` with `p,q ∈ {s,c,d,n}` are also
provided, but this function is more efficient if more than one elliptic
function with the same arguments is required.
"""
function ellipj end


#-----------------------
# Convenience functions

chars = ("s","c","d")
for (i,p) in enumerate(chars)
    pn = Symbol(p*"n")
    np = Symbol("n"*p)
    @eval begin
        $pn(u,m) = ellipj(u,m)[$i]
        $np(u,m) = 1/$pn(u,m)
    end
end
for p in (chars...,"n")
    pp = Symbol(p*p)
    @eval $pp(u,m) = one(promote_type(typeof.(float.((u,m)))...))
end

for p in chars, q in chars
    p == q && continue
    pq = Symbol(p*q)
    pn = Symbol(p*"n")
    qn = Symbol(q*"n")
    @eval begin
        function $pq(u::Number,m::Number)
            sn,cn,dn = ellipj(u,m)
            return $pn/$qn
        end
    end
end

for p in (chars...,"n"), q in (chars...,"n")
    pq = Symbol(p*q)
    @eval begin
"""
    $(string($pq))(u,m)
Jacobi elliptic function `$($p)$($q)`.
See also `ellipj(u,m)` if more than one Jacobi elliptic function
with the same arguments is required.
"""
function $pq end
    end
end

"""
    am(u,m)
Jacobi amplitude function
"""
function am(u, m)
    asin(sn(u, m))
end


end #Jacobi module
