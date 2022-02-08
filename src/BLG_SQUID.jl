module BLG_SQUID
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) BLG_SQUID

    using Requires

    function __init__()
        @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" begin       
        @require  VegaLite = "112f6efa-9a02-5b7d-90c0-432ed331239a" begin 
        end
        end
    end

    using SharedArrays, Distributed
    using Quantica, StaticArrays, Parameters, LinearAlgebra, StatsBase
    using Baselet, Arpack, CSV, DataFrames, Dates
    using Optim, ProgressMeter, Interpolations, Random, LsqFit
    
    using PhysicalConstants.CODATA2018: ustrip, @u_str, h, ħ, k_B, m_e, e, μ_B
    
    export Params
    export nanoribbonS, nanoribbonSA, nanoribbonSZ, Params,  modelS, rectangle_weaklink, 
        rectangle_randombounds_sc, ldosonlattice_averaged_sc, ldosonlattice
    export icϕ_exactdiag, fraunhofer_abs_exact, maxicϕ_exactdiag, fraunhofer_abs_exact_adaptive

    include("model.jl")
    include("nanoribbon.jl")
    include("bounded_sys.jl")
    include("ldos.jl")
    include("spectrum.jl")
    include("weaklink.jl")
    include("save.jl")
    include("adaptive_maxfinder.jl")
    include("transport.jl")

end
