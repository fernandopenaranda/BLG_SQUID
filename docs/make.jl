using BLG_SQUID
using Documenter

DocMeta.setdocmeta!(BLG_SQUID, :DocTestSetup, :(using BLG_SQUID); recursive=true)

makedocs(;
    modules=[BLG_SQUID],
    authors="Fernando Peñaranda",
    repo="https://github.com/fernandopenaranda/BLG_SQUID.jl/blob/{commit}{path}#{line}",
    sitename="BLG_SQUID.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fernandopenaranda.github.io/BLG_SQUID.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/BLG_SQUID.jl",
    devbranch="main",
)
