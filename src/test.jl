using Plots
"""
Datoteka s kodo za določanje št. razdelitev intervala za računanje vrednosti integrala funkcije f(x) = sin(x)/x
"""

"""
test_kvadrature()

Testira racunanje integrala funkcije f(x) = sin(x)/x z uporabo Gauss-Lengendrovih kvadratur in ugotavlja, na koliko
delov `n` je potrebno interval razdeliti, da napaka pade pod eps=1e-10.
Za natančno vrednost integrala uporabi vrednost, ki jo izračuna Wolfram Alpha.

# Primer
```jldoctest
julia> test_kvadrature()
``` 
"""
function test_kvadrature()
    # natancna resitev iz Wolfram Alpha
    resitev = 1.549931244944674137274408400730639012183184893966372210477969710681487208951511074986007223927691325
    
    eps = 1e-10
    err = 1
    errors = []

    f(x) = sin.(x)./x
    a = eps
    b = 5.0

    n=10
    razdelitve = []
    while err > eps && n < 1000000000

        I = gauss_legendre(f, n, a, b)

        err = abs(I - resitev)
        push!(errors, err)
        push!(razdelitve, n)
        n += 10
        println("N: ", n, "; err: ", err)
    end

    p = plot(razdelitve, errors)
    savefig(p, "integracija_napake.png")

    return razdelitve, errors
end


```
julia> test_kvadrature()
N: 20; err: 2.1702422339231475e-6
N: 30; err: 1.348454383709452e-7
N: 40; err: 2.6530374563904502e-8
N: 50; err: 8.323194133907919e-9
N: 60; err: 3.349605925961896e-9
N: 70; err: 1.5634429306743414e-9
N: 80; err: 7.97840238320191e-10
N: 90; err: 4.262803443566554e-10
N: 100; err: 2.2854607095723622e-10
N: 110; err: 1.1555578716127002e-10
N: 120; err: 4.7226889066109834e-11

```