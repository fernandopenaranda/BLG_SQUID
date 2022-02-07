"computes the fraunhofer pattern only considering the states inside the parent gap using 
exact diagonalization"
fraunhofer_abs_exact(Bperp::Array{T,1}, p; kw...) where {T} = 
    fraunhofer_abs_exact(Bperp, p, missing; kw...) 
function fraunhofer_abs_exact(Bperp::Array{T,1}, p,
        Δϕ::Union{Array{T,1}, Missing}; kw...) where {T}
    println("computing current matrix...")
    icmax = SharedArray(similar(Bperp))
    ic = icϕ_exactdiag(Bperp, p, Δϕ; kw...)
    for i in 1:length(Bperp)
        println("it: ", i/length(Bperp))
        icmax[i] = maximum(abs.(ic[:, i]))
    end
    return Bperp, icmax, ic
end


function fraunhofer_abs_exact_adaptive(Bperplist, p; kw...)
    icmax = maxicϕ_exactdiag(Bperplist, p; kw...) 
    return Bperp, icmax
end

"""
    `icϕ_exactdiag(Bperplist::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing}; kw...)`
supercurrent sweep with both Bperp and Δϕ. See: `supercurrent_exactdiag()`
"""
function icϕ_exactdiag(Bperplist::Array{T,1}, p, Δϕlist; kw...) where {T}
    # Δϕlist = ifelse(isa(Δϕ, Missing) == true, -0.5:1:π+0.5, Δϕ)
    I = zeros(Float64, length(Δϕlist), length(Bperplist))
    [I[:, i] = supercurrent_exactdiag(collect(Δϕlist), 
        reconstruct(p, B = SA[0,0,Bperplist[i]]); kw...) for i in 1:length(Bperplist)]
    return I 
end

"""
    `maxicϕ_exactdiag(Bperplist::Array{T,1}, p, Δϕ::Union{Array{T,1}, Missing}; kw...)`
computes the Ic_optim = max(Ic, ϕ) using a ML K-cross validation search algorithm for the 
minimum. See: `adaptive_max_finder.jl`
"""
function maxicϕ_exactdiag(Bperplist::Array{T,1}, p; kw...) where {T}
    icmax = SharedArray(similar(Bperplist))
    @sync @distributed for i in 1:length(Bperplist) 
        hpar = 
            rectangle_weaklink(reconstruct(p, B = SA[0,0, Bperplist[i]]), false)

        gen_model(x::Array, hp; kw...) = 
            supercurrent_exactdiag_adaptive(x, hp; kw...)  #not sure if p is being updated as it should
        gen_model(x::Array, f::Array, xmax::Number, hp; kw...) = 
            supercurrent_exactdiag_adaptive(x, f, xmax, hp; kw...)
        
        icmax[i] = adaptive_max_finder(gen_model, fourier_model, hpar; kw)[1]
    end
    return icmax 
end

"""
    `supercurrent_exactdiag(Δϕlist, p = Params() ; nev = 10)` 
Computes the supercurrent for a fixed Bperp value specified in p. If `nev::Missing` it 
iteratively finds all ABS (i.e. with energies inside a lengthscale which we set us the gap).
It returns the supercurrent contribution of these subgap subspace. Note that this subspace
is chosen to always capture all the ABS contribution although it may contain also continuum
states, because the parent gap will be larger that the induced gap. However, once the
subspace is specified, we can always substract this contribution for the full KPM 
contribution. 
    `method = :free_energy`
"""
supercurrent_exactdiag(Δϕlist, p = Params()::Params; kw...) = 
    supercurrent_exactdiag(Δϕlist, rectangle_weaklink(p, false); kw...)

function supercurrent_exactdiag(Δϕlist, hpar::Quantica.ParametricHamiltonian; kw...)
    f = SharedArray(zeros(Float64, length(Δϕlist)))
    @sync @distributed for i in 1:length(Δϕlist)
        println(i/length(Δϕlist))
        f[i] = f_e(hpar, Δϕlist[i]; kw...)
    end
    ip = interpolate((Δϕlist,), f, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in Δϕlist]
    return @. 2/ħoec .* deriv
end

"""
`supercurrent_exactdiag_adaptive()`
Same as `supercurrent_exactdiag()` but for the fact it uses the adaptive_max_finder approach
instead of a sweep in SC phase differences.
 See: `adaptive_max_finder.jl`
"""

function supercurrent_exactdiag_adaptive(ϕ, f, ϕnew, hpar::Quantica.ParametricHamiltonian; kw...)
    fnew = f_e(hpar, ϕnew; kw...)
    ϕn = vcat(ϕ, ϕnew)
    fn = vcat(f, fnew)
    perm = sortperm(ϕn)
    ϕn = ϕn[perm]
    fn = fn[perm]
    newindex = findall(x-> x == ϕnew, ϕn)[1]
    ip = interpolate((ϕn,), fn, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in ϕn] # check if the returned index Ic# is the correct one
    return fnew, @. 2/ħoec * deriv[newindex]
end    

function supercurrent_exactdiag_adaptive(Δϕlist, hpar::Quantica.ParametricHamiltonian; kw...)
    f = zeros(Float64, length(Δϕlist))
    for i in 1:length(Δϕlist)
        println("init phases: ", i/length(Δϕlist))
        f[i] = f_e(hpar, Δϕlist[i]; kw...)
    end
    ip = interpolate((Δϕlist,), f, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in Δϕlist]
    return Δϕlist, f, @. 2/ħoec .* deriv
end

"""
returns the free energy computed using the exact eigenvalues for the exact evaluation of
the supercurrent over a finite set of eigenvectors
""" 
f_e(hpar::Quantica.ParametricHamiltonian, Δϕ; kw...) = -sum(negative_eigen(hpar, Δϕ))
    
"""
returns the NEGATIVE eigenvalues withing the gap using for the exact supercurrent
calculation an adaptive procedure for setting the number of nevs calculated with the shift
invert method
"""
function negative_eigen(hpar, Δϕ)
    sp = spectrum(hpar(ϕ = Δϕ), method = ArpackPackage(nev = 104, sigma = -0.001im))
    # sp = DACPdiagonaliser(hpar(ϕ = Δϕ), 0.1)
    λ = real(sp.energies)
    λneg = λ[λ.<=0];
    return λneg
end