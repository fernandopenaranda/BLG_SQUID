using StaticArrays

p = Params(Ln = 200, W = 200, Ls = 20, Ws = 20, scale = 40, λ = 5,
          Δ = 0, d = 0, τ = 1, B = SA[0,0,.0], EZ = SA[0,0,0], μN = 1e-6);
ph = rectangle_weaklink(p, false);

fΦ = 1
Blist = collect(0.:fΦ/2:5*fΦ)
Δϕ =  -0.5:3:π+0.5
deltaphi = collect(Δϕ).*0

icϕ_exactdiag(Blist, p, deltaphi)
# vlplot(ph(ϕ=0), maxdiameter = 4, plotlinks = false)
maxicϕ_exactdiag(Blist, p)

# p = Params(Ln = 195, W = 2500, Ls = 20, Ws = 20, scale = 40, λ = 5,
#           Δ = 1, d = 0, τ = 1, B = SA[0,0,1], EZ = SA[0,0,0], μN = 1e-6);



#####

p1 = Params(Ln =50, W = 50, Ls = 3, Ws = 3, scale = 40, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,.0], EZ = SA[0,0,0], μN = 1e-6);
# ph1 = rectangle_weaklink(p1, false);


# p1 = Params(Ln = 200, W = 200, Ls = 3, Ws = 3, scale = 40, λ = 5,
#           Δ = 1, d = 0, τ = 1, B = SA[0,0,.0], EZ = SA[0,0,0], μN = 1e-6);
# ph1 = rectangle_weaklink(p1, false);

# p2 = Params(Ln = 200, W = 200, Ls = 3, Ws = 3, scale = 40, λ = 5, U = 10,
#           Δ = 1, d = 0, τ = 1, B = SA[0,0,.0], EZ = SA[0,0,0], μN = 1e-6);
# ph2 = rectangle_weaklink(p2, false);

# sp1 = spectrum(ph1(ϕ = 1),  method = ArpackPackage(nev = 52, sigma = -0.001im))
# en1 = sort(abs.(sp1.energies))

# sp2 = spectrum(ph2(ϕ = 1),  method = ArpackPackage(nev = 52, sigma = -0.001im))
# en2 = sort(abs.(sp2.energies))

