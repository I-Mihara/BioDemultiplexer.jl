using BioDemultiplexer
using Documenter

DocMeta.setdocmeta!(BioDemultiplexer, :DocTestSetup, :(using BioDemultiplexer); recursive=true)

makedocs(;
    modules=[BioDemultiplexer],
    authors="I.Mihara <issei.mihara@xforestx.com> and contributors",
    repo="https://github.com/I-Mihara/BioDemultiplexer.jl/blob/{commit}{path}#{line}",
    sitename="BioDemultiplexer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://I-Mihara.github.io/BioDemultiplexer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/I-Mihara/BioDemultiplexer.jl",
    devbranch="main",
)
