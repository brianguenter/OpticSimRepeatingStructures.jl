# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Documenter
using OpticSimRepeatingStructures
import Luxor
import OpticSim

# override certain functions to allow production of interactive figures
OpticSim.set_current_mode(:docs)


makedocs(
    sitename="OpticSimRepeatingStructures.jl",
    format=Documenter.HTML(
        # prettyurls = get(ENV, "CI", nothing) == "true",
        assets=[asset("assets/logo.svg", class=:ico, islocal=true)],
    ),
    modules=[OpticSim],
    pages=[
        "Home" => "repeat.md"
    ]
)

deploydocs(
    repo="github.com/brianguenter/OpticSimRepeatingStructures.jl.git",
    devbranch="main",
    push_preview=true,
)

