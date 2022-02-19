# p = Params(Ln = 400, W = 400, Ls = 40, Ws = 20, scale = 40, λ = 5., U = 0, α = 0,
#     Δ = 1, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,1e-6], μN = 1e-6);


# p = Params(Ln = 50, W = 80, Ls = 10, Ws = 0, scale = 20, λ = 5., U = 0, α = 0,
# Δ = 5, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,1e-6], μN = 1e-6);

p = Params(Ln = 200, W = 320, Ls = 30, Ws = 0, scale = 20, λ = 5., U = 0, α = 0,
Δ = 5, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,1e-6], μN = 1e-6);

fΦ = 2π *BLG_SQUID.quantumflux(p)
Blist = collect(0.:fΦ/4:5*fΦ)
 
# ph = rectangle_weaklink(p, false, true, Δx_mask = 20, Δy_mask = 20)
# ic = maxicϕ_exactdiag(Blist, p, mask=true, Δx_mask = 30, Δy_mask = 30)
ic = icϕ_exactdiag(Blist, p, mask=true, Δx_mask = 20, Δy_mask = 20)

@time sp = spectrum(ph(ϕ = 0), method = ArpackPackage(nev = 204, sigma = 0.0im))

function ldosonlattice(number) 
    println(sp.energies[number])
   return ldosonlattice(sp.states[:,number], ph(ϕ = 0))
end



ph = rectangle_weaklink(p, false, false)

@time sp = spectrum(ph(ϕ = 0), method = ArpackPackage(nev = 192, sigma = 0.0im))



##########
# check regime

function Bcriteria(B::Array, scale; threshold = 50.)
    ħc = ustrip(u"meV*nm",ħ*c_0)
    a0 = 0.246    #nm
    Bmax = findmax(abs.(B))[1]
    l_B = ħc/Bmax
    l_scale = scale * a0    
    println(l_B/ l_scale)
    println(l_B/a0)
    l_B/ l_scale > threshold ? nothing : @warn("lB ≈ a0*scale, decrease scaling")
    l_B/ a0 > threshold ? nothing : @warn("lB ≈ a0, increase BLG area (we are in t
        he Hofstadter regime)")
end

#######################
# 
p = Params(Ln = 10, W = 10, Ls = 10, Ws = 0, scale = 10, λ = 5., U = 0, α = 0,
Δ = 5, d = 0, τ = 1, B = SA[0,0,0], EZ = SA[0,0,1e-6], μN = 1e-6);
Blist = [0.]
ic = icϕ_exactdiag(Blist, p, mask = true, Δx_mask = 30, Δy_mask = 30)