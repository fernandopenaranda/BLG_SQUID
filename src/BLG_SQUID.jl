module BLG_SQUID
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) BLG_SQUID

    using Requires

 
    using CairoMakie, VegaLite 
    using SharedArrays, Distributed
    using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
    using Baselet, Arpack, CSV, DataFrames, Dates
    using Optim, ProgressMeter, Interpolations, Random, LsqFit
    
    using PhysicalConstants.CODATA2018: ustrip, @u_str, h, ħ, k_B, m_e, e, μ_B
    
    export Params
    export nanoribbonS, nanoribbonSA, nanoribbonSZ, Params,  modelS, rectangle_weaklink, 
        rectangle_randombounds_sc, ldosonlattice_averaged_sc, ldosonlattice
    export icϕ_exactdiag, fraunhofer_abs_exact, maxicϕ_exactdiag, fraunhofer_abs_exact_adaptive

    export Params_eff_quad, supercurrent_ef, supercurrent2_1dchains, spectrumvsphase_ef,
        fraunhofer_eff, plotfraunhof_vs_phase 

    include("model.jl")
    include("nanoribbon.jl")
    include("bounded_sys.jl")
    include("ldos.jl")
    include("spectrum.jl")
    include("weaklink.jl")
    include("save.jl")
    include("adaptive_maxfinder.jl")
    include("transport.jl")
    include("effective_model.jl")

end
