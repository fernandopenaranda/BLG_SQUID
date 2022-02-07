using StaticArrays

p = Params(Ln = 100, W = 1000, Ls = 3, Ws = 3, scale = 40, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,.0], EZ = SA[0,0,0], μN = 1e-6);
ph = rectangle_weaklink(p, false);

Blist = [1e-5,0.05]
Δϕ =  -0.5:1:π+0.5

icϕ_exactdiag(Blist, p, collect(Δϕ))
# vlplot(ph(ϕ=0), maxdiameter = 4, plotlinks = false)
maxicϕ_exactdiag(Blist, p)

# p = Params(Ln = 195, W = 2500, Ls = 20, Ws = 20, scale = 40, λ = 5,
#           Δ = 1, d = 0, τ = 1, B = SA[0,0,1], EZ = SA[0,0,0], μN = 1e-6);
