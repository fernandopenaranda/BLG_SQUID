using StaticArrays

p = Params(Ln = 200, W = 200, Ls = 20, Ws = 20, scale = 40, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,fΦ], EZ = SA[0,0,0], μN = 1e-6);
fΦ = 0.05
Blist = collect(0.:fΦ/6:10*fΦ)
ic = maxicϕ_exactdiag(Blist, p)
# ph = rectangle_weaklink(p, false);
##############################
icϕ_exactdiag(Blist, p, deltaphi)

#####

p = Params(Ln =50, W = 50, Ls = 20, Ws = 20, scale = 40, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,0.025], EZ = SA[0,0, 0], μN = 1e-6);
ph = rectangle_weaklink(p, false);
shortlist = [0., 0.005]
ic = maxicϕ_exactdiag(shortlist, p)
