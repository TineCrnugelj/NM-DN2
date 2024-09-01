module DN2

"""Koda za približno računanje integrala, ki uporablja Gauss-Legendrove kvadrature
"""

export DN2


"""
I = gauss_legendre(f, n, a, b)

Izračuna določeni integral funkcije `f` na intervalu od `[a, b]`, pri čemer
interval razdeli na `n` delov.
Za računanje integrala uporablja Gauss-Legendrove kvadrature na 2 točkah.

# Primer
```jldoctest
julia> f(x) = sin.(x)./x
julia> I = gauss_legendre(f, 100, 0.0, 5.0)
1.5499312451602307
``` 
"""
function gauss_legendre(fn::Function, n::Int, a::Float64, b::Float64)::Float64
    w_i = [1, 1]    # uteži
    x_i = [-1/sqrt(3), 1/sqrt(3)]   # x vrednosti

    I = 0.0
    h = (b - a)/n   # del intervala
    
    for i=1:n
        I += w_i' * fn(h/2 .* x_i .+ (a + h*(i - 1/2)))
    end

    return h/2 * I
end

end