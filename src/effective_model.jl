const τ0 = @SMatrix[1 0; 0 1]
const τx = @SMatrix[0 1; 1 0]
const τy = @SMatrix[0 -im; im 0]
const τz = @SMatrix[1 0; 0 -1]

const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σ0τ0 = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σzτ0 = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
const σzτz = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
const σyτy = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]
const σyτz = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]
const σyτ0 = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 -im; 0 0 im 0]
const σxτz = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]
const σxτ0 = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]

@with_kw struct Params_eff_quad @deftype Float64
    μN = 0
    Δ = 1
    t0 = 1
    a0 = 1
    Ln = 5
    Ls = 10
    psweight = 1/20
end

@with_kw struct Params_eff_linear @deftype Float64
    μN = 0
    Δ = 1
    t0 = 4
    a0 = 1
    Ln = 5
    Ls = 10
    psweight = 1/20
end

function h_eff_cuadratic(p = Params_eff_quad())
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(-φ)])

    lat =  unitcell(LatticePresets.linear(a0 = a0), region = r -> abs(r[1]) < Ln + Ls)
    
    mod_onsite = onsite(- μN * τz)
    mod_hop = hopping(-t0 * τz, range = a0)
    mod_on! = @onsite!((o; θ) -> o + θ*psweight * τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * diagphi(sign(dr[1])*θ))
    self_region(r) = !(-Ln <= r[1] <= Ln)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o + self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+θ)/4) * τy * diagphi(sign(r[1])*(ϕ+θ)/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end


function h_eff_cuadratic_unbounded(p = Params_eff_quad())
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(-φ)])

    lat =  LatticePresets.linear(a0 = a0)
    
    mod_onsite = onsite(- μN * τz)
    mod_hop = hopping(-t0 * τz, range = a0)
    mod_on! = @onsite!((o; θ) -> o + θ*psweight * τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * diagphi(sign(dr[1])*θ))
    self_region(r) = !(abs.(r[1]) < 0)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o + self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+θ)/4) * τy * diagphi(sign(r[1])*(ϕ+θ)/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end

function h_eff_helical(p = Params_eff_linear())
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    lat =  unitcell(LatticePresets.linear(a0 = a0), region = r -> abs(r[1]) < Ln + Ls)
    
    mod_onsite = onsite(0I)
    mod_hop = hopping((r, dr) -> -1im * t0 * sign(dr[1]) * σ0τz, range = a0)
    mod_on! = @onsite!((o, r; θ) -> o + θ * psweight * σ0τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * 1 + 0 * diagphi(sign(dr[1])*θ))
    self_region(r) = !(-Ln <= r[1] <= Ln)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o + self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+θ)/4) * σyτy * diagphi(sign(r[1])*(ϕ+θ)/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown,:hup, :hdown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end

function h_eff_helical_unbounded(p = Params_eff_linear())
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    lat =  LatticePresets.linear(a0 = a0)
    
    mod_onsite = onsite(0I)
    mod_hop = hopping((r, dr) -> -1im*t0 * sign(dr[1]) * σ0τz, range = a0)
    mod_on! = @onsite!((o; θ) -> o + θ * psweight * σ0τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * 1 + 0 * diagphi(sign(dr[1])*θ))
    self_region(r) =  !(abs.(r[1]) < 0)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o +self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+θ)/4) * σyτy * diagphi(sign(r[1])*(ϕ+θ)/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown,:hup, :hdown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end

function h_eff_helical_coupled(p = Params_eff_linear(Ln = 0.5, Ls =1), )
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    lat =  unitcell(LatticePresets.square(a0 = a0), 
        region = r -> (abs(r[1]) < Ln + Ls &&  0 <= r[2] <= 3a0 ))
    ymin, ymax = extrema(r->r[2], sitepositions(lat))
    mod_onsite = onsite(0I)
    function hop(r, dr)
        if dr[2] == 0.
            -1im * sign(dr[1])
        elseif dr[2] != 0 && abs(r[1]) <= Ln -a0/2
            0
        else 
            0.1
        end
    end
    mod_hop = hopping((r, dr) -> hop(r, dr) * t0  * σ0τz, range = a0)
    mod_on! = @onsite!((o, r; θ) -> o + θ * r[2]/abs(ymax-ymin) * psweight * σ0τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * 1 + 0 * diagphi(sign(dr[1])*θ))
    self_region(r) = !(-Ln < r[1] < Ln)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o + self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+(θ*r[2]/abs(ymax-ymin)))/4) * σyτy * diagphi(sign(r[1])*(ϕ+(θ*r[2]/abs(ymax-ymin)))/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown,:hup, :hdown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end

function h_eff_helical_coupled_gaugex(p = Params_eff_linear(Ln = 0.5, Ls =1))
    (; μN, Δ, t0, a0, Ln, Ls, psweight) = p
    diagphi(φ) = Diagonal(SA[cis(φ), cis(φ), cis(-φ), cis(-φ)])
    lat =  unitcell(LatticePresets.square(a0 = a0), 
        region = r -> (abs(r[1]) < Ln + Ls &&  0 <= r[2] <= 1a0 ))
    mod_onsite = onsite(0I)
    function hop(r, dr)
        if dr[2] == 0.
            -1im * sign(dr[1])
        elseif dr[2] != 0 && abs(r[1]) <= Ln -a0/2
            0
        else 
            0
        end
    end
    mod_hop = hopping((r, dr) -> hop(r, dr) * t0  * σ0τz, range = a0)
    mod_on! = @onsite!((o, r; θ) -> o + θ * r[1] * psweight * σ0τz) #ps need to be renormalized
    mod_peierls! = @hopping!((t, r, dr; θ) -> t * 1 + 0 * diagphi(sign(dr[1])*θ))
    self_region(r) = !(-Ln < r[1] < Ln)
    mod_sc! = @onsite!((o, r; ϕ, θ) -> o + self_region(r) * Δ * 
        diagphi(sign(r[1])* (ϕ+θ)/4) * σyτy * diagphi(sign(r[1])*(ϕ+θ)/4)') 

    ham = lat |> hamiltonian(mod_onsite + mod_hop, orbitals = (:eup, :edown,:hup, :hdown)) |>
         parametric(mod_on!, mod_sc!, mod_peierls!)
    return ham
end

function spectrumvsphase_ef(p = Params_eff_quad(), Δϕlist = collect(-π:π/100:π), θ = 0)
    hp = h_eff_helical_coupled(p)#h_eff_helical(p)#h_eff_cuadratic(p)
    s = spectrum(hp(ϕ = 0, θ = 0)).energies
    sp = zeros(Float64, length(s), length(Δϕlist))
    for i in 1:length(Δϕlist)
        sp[:, i] = real( spectrum(hp(ϕ = Δϕlist[i], θ = θ)).energies)
    end
    return sp
end

supercurrent_ef(p::Params_eff_quad, Δϕlist = collect(-2π:π/100:2π), θ = 0) = 
    supercurrent_ef(spectrumvsphase_ef(p, Δϕlist, θ), Δϕlist)

function supercurrent_ef(sp::Matrix{Float64}, Δϕlist = collect(-2π:π/100:2π))
    f = fe(sp)
    ip = interpolate((Δϕlist,), f, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in Δϕlist]
    return @. 2/ħoec .* deriv
end

fe(sp::Matrix) = [fe(sp[:,i]) for i in 1:size(sp,2)]
fe(sp::Vector) = -sum(sp[sp.≤0])


function supercurrent2_1dchains(θ, p = Params_eff_quad(), Δϕlist = collect(-π:π/100:π))
    i1 = supercurrent_ef(spectrumvsphase_ef(p, Δϕlist, θ), Δϕlist)
    # i2 = supercurrent_ef(spectrumvsphase_ef(p, Δϕlist, -θ), Δϕlist)
    # ploti1i2_vs_phase(i1,i2, Δϕlist)
    return i1#+i2
end

function fraunhofer_eff(p = Params_eff_quad(), θlist = 0:π/3:π, Δϕlist = collect(-2π:π/100:2π))
    pattern = [maximum(abs.(supercurrent2_1dchains(θlist[i], p, Δϕlist))) for i in 1:length(θlist)]
    plotfraunhof_vs_phase(pattern, θlist)
    # return pattern
end

###
using CairoMakie
function plotsp_vs_phase(sp, Δϕ, ylims = (-1.2,1.2))
    f = Figure(resolution = (500, 500))
    ax = Axis(f[1,1], xlabel= "ϕ/π", ylabel = "E (meV)")
    for i in 1:size(sp,1)
        scatter!(ax, Δϕ/π, sp[i, :], color = :gray, markersize = 2)
        lines!(ax, Δϕ/π, sp[i, :], color = :gray)
    end
    return f
end

function ploti_vs_phase(I, Δϕ)
    f = Figure(resolution = (500, 500))
    ax = Axis(f[1,1], xlabel= "ϕ/π", ylabel = "I (arb. units)")
    lines!(ax, Δϕ/π, I)
    return f    
end

function ploti1i2_vs_phase(i1,i2, Δϕ)
    f = Figure(resolution = (500, 500))
    ax = Axis(f[1,1], xlabel= "ϕ/π", ylabel = "I (arb. units)")
    lines!(ax, Δϕ/π, i1, color = :red)
    lines!(ax, Δϕ/π, i2, color = :blue)
    return f    
end

function plotfraunhof_vs_phase(pattern, θlist)
    f = Figure(resolution = (500, 500))
    ax = Axis(f[1,1], xlabel= "θ", ylabel = "Ic (arb. units)")
    lines!(ax, θlist/π, pattern)
    return f    
end

###

# user
# es = spectrumvsphase_ef();
# ic = supercurrent_ef();
# plotsp_vs_phase(es, collect(-π:π/100:π))
# ploti_vs_phase(ic, collect(-2π:π/100:2π))
# BLG_SQUID.isnambu_spinless( BLG_SQUID.h_eff_cuadratic())

function isnambu_spinless(h::Quantica.ParametricHamiltonian) 
    isnambu_spinless(h(ϕ = 1, θ = 1))
end
function isnambu_spinless(h::Quantica.Hamiltonian)
    
    ham =  Quantica.flatten(h).harmonics[1].h 
    dim = size(ham,1)
    dim2 = Int64(size(ham,1)/2)
    basis_change_mat = zeros(ComplexF64, dim, dim)
    step = 2
    array_indices = sort(collect(range(1, stop = Int64(dim), step = step)))
    [basis_change_mat[i, array_indices[i]] = 1. for i in 1:1:length(array_indices)]
    array_indices_h = sort(collect(range(2, stop = Int64(dim), step = step)))
    [basis_change_mat[Int64(dim/2)+i, array_indices_h[i]] = 1.
        for i in 1:1:length(array_indices_h)]
    bs = basis_change_mat * Quantica.flatten(h).harmonics[1].h * basis_change_mat'
    l =  -conj.(bs[1:dim2,1:dim2]) == bs[dim2+1:end,dim2+1:end]
    if l == true
        return l
    else
        println("false")
    return -conj.(bs[1:dim2,1:dim2]) - bs[dim2+1:end,dim2+1:end]
    end
end