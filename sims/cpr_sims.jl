using StaticArrays

p = Params(Ln = 200, W = 200, Ls = 20, Ws = 20, scale = 40, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,0], μN = 1e-6);
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


p = Params(Ln = 200, W = 500, Ls = 20, Ws = 20, scale = 20, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,0], μN = 1e-6);

BLG_SQUID.quantumflux(p)
fΦ = 0.02067
Blist = collect(0.:fΦ/2:3*fΦ)
ic = maxicϕ_exactdiag(Blist, p, mask=true, Δx_mask = 40, Δy_mask = 40)
# ph = rectangle_weaklink(p, false, true, Δx_mask = 20, Δy_mask = 20)



p = Params(Ln = 50, W = 80, Ls = 10, Ws = 10, scale = 10, λ = 5,
          Δ = 1, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,0], μN = 1e-6);

fΦ = BLG_SQUID.quantumflux(p)
Blist = [0.]#collect(0.:fΦ/2:3*fΦ)
ic = maxicϕ_exactdiag(Blist, p, mask=true, Δx_mask = 1000, Δy_mask = 1000)


h = BLG_SQUID.nanoribbonN(p, (0,1));
b = bandstructure(h, 
            cuboid((-.5, .5); subticks = 151), 
            mapping = x -> 2pi*x, 
            method = ArpackPackage(nev=32, maxiter = 300, sigma=-0.0000));