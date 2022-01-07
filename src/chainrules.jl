# Rules for Elliptic E
using ChainRulesCore

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

# function ChainRulesCore.rrule(::typeof(Jacobi.am), u, m)
#     y = am(u, m)
#     function am_pullback(ȳ)
#         f̄ = NoTangent()
#         ū = dn(u, m)
#         m̄ = @thunk begin
#             central
#         end
#     end
# end
