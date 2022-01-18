# Rules for Elliptic E
using ChainRulesCore
using ForwardDiff

function ChainRulesCore.rrule(::typeof(Elliptic.E), phi, m)
    y = E(phi, m)
    function E_pullback(ȳ)
        f̄ = NoTangent()
        p̄hi = @thunk(sqrt(1-m^2*sin(phi)^2)*ȳ)
        m̄ = @thunk(1/(2*m)*(y - F(phi, m))*ȳ)
        return f̄, p̄hi, m̄
    end
    return y, E_pullback
end

function ChainRulesCore.rrule(::typeof(Elliptic.F), phi, m)
    y = F(phi, m)
    function F_pullback(ȳ)
        f̄ = NoTangent()
        p̄hi = @thunk(ȳ/sqrt(1-m^2*sin(phi)^2))
        m̄ = @thunk begin
                t1 = E(phi, m)/(2*(1-m)*m)
                t2 = -y/(2*m)
                t3 = -sin(2*phi)/(4*(m-1)sqrt(1-m*sin(phi)^2))
                (t1+t2+t3)*ȳ
        end
        return f̄, p̄hi, m̄
    end
    return y, F_pullback
end

function ChainRulesCore.rrule(::typeof(K), m)
    k,e = ellipke(m)
    function ellipke_pullback(ȳ)
        f̄ = NoTangent()
        m̄ = @thunk((e - (1-m)*k)/(2*(1-m)*m)*ȳ)
        return f̄, m̄
    end
    return k, ellipke_pullback
end

function ChainRulesCore.rrule(::typeof(E), m)
    k,e = ellipke(m)
    function ellipke_pullback(ȳ)
        f̄ = NoTangent()
        m̄ = @thunk((e-k)/(2m)*ȳ)
        return f̄, m̄
    end
    return k, ellipke_pullback
end

function ChainRulesCore.rrule(::typeof(Jacobi.sn), u, m)
    sn, cn, dn = ellipj(u, m)
    ms(x) = Jacobi.sn(u, x)
    function sn_pullback(ȳ)
        f̄ = NoTangent()
        ū = cn*dn
        m̄ = ForwardDiff.derivative(ms, m)
        return  f̄, ū, m̄
    end
    return sn, sn_pullback
end

function ChainRulesCore.rrule(::typeof(Jacobi.cn), u, m)
    sn, cn, dn = ellipj(u, m)
    ms(x) = Jacobi.cn(u, x)
    function cn_pullback(ȳ)
        f̄ = NoTangent()
        ū = -sn*dn
        m̄ = ForwardDiff.derivative(ms, m)
        return  f̄, ū, m̄
    end
    return cn, cn_pullback
end

function ChainRulesCore.rrule(::typeof(Jacobi.dn), u, m)
    sn, cn, dn = ellipj(u, m)
    ms(x) = Jacobi.dn(u, x)
    function dn_pullback(ȳ)
        f̄ = NoTangent()
        ū = -m*sn*cn
        m̄ = ForwardDiff.derivative(ms, m)
        return  f̄, ū, m̄
    end
    return dn, dn_pullback
end
