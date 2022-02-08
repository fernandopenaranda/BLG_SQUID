"""
Implements configuration B, i.e. BLG in Bernal stacking with self energy on both the ZZ and
the AC edges. There is a different phase difference between left and right sides of the 
weak link (-ϕ/2 for x<0 ϕ/2 for x>0).
Neither smoothness nor leads are implemented. Also disorder and rotation are not considered.
"""
#mask == false ?  devicereg = reg : devicereg = xor(reg, 
    #RegionPresets.rectangle((Ln-2*Δx_mask, W-2*Δy_mask), (0., W/2 + a0 *scale)))
function rectangle_weaklink(p, selfy = true)
    (; Ls, Δ, Ws, B, Ln, W) = p
    (; model0, field!, modelinter) = modelS(p)
    lat_top, lat_bot = latBLG(p, 0) #0 angle rotation
    
    println("flux: ",B[3]*Ln*W/Φ0)
    # PEIERLS PHASES
    diagphi(φ) = Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    piecewise(y) = clamp(y, -W/2, W/2) # clamp(x, -Ln/2 - a0/(2*√3) , Ln/2 + a0/(2*√3)))
    A(r, B) = [-B[3]/ħoec, 0, 0] * piecewise(r[2])
    eφ(r, dr, B) = diagphi(dot(A(r, B), dr))
    peierls! =@hopping!((t, r, dr) -> t * eφ(r, dr, B))
    
    # LOCAL SELF-ENERGY MODEL
    xmin, xmax = extrema(r->r[1], sitepositions(Quantica.combine(lat_top, lat_bot)))
    ymin, ymax = extrema(r->r[2], sitepositions(Quantica.combine(lat_top, lat_bot)))
    self_region(r) = ifelse(selfy == false, !(xmin+Ls <= r[1] <= xmax - Ls),  
        !(xmin+Ls <= r[1] <= xmax - Ls) || !(ymin+Ws <= r[2] <= ymax - Ws))

    sCself! = @onsite!((o, r; ϕ) -> o + self_region(r) * Δ *
        eφ(r, r+[0,0,0], B) * diagphi(sign(r[1])* ϕ/4) * σyτy *
            diagphi(sign(r[1])*ϕ/4)' * eφ(r, r+[0,0,0], B)')
              
    # HAMILTONIAN BUILD
    h_top = lat_top |> hamiltonian(model0; orbitals = Val(4)) |> 
        unitcell(mincoordination = 5)    
    h_bot = lat_bot |> hamiltonian(model0; orbitals = Val(4)) |>
        unitcell(mincoordination = 5)
    ph = Quantica.combine(h_top, h_bot; coupling = modelinter) |> 
        parametric(peierls!, sCself!) 
            return ph
end

"""
    spectrumvsphase(philist, p; nev = 16, kw...)
computes a number nev of eigenvalues as a function of the phase difference across the weak
link. 
    see: rectangle_weakling(::Params()), Params(), modelS()
"""
function spectrumvsphase(philist, p; nev = 16, kw...)
    ph = rectangle_weaklink(p; kw...)
    elist = zeros(Float64, nev, length(philist))
    for i in 1:length(philist)
        s = spectrum(ph(phi = philist[i]), method = ArpackPackage(sigma = 1e-7, nev = nev))
        elist[:,i] = s.energies
    end
    return philist, elist
end
