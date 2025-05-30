# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Documenter
using OpticSimRepeatingStructures
import Luxor
import OpticSim




makedocs(
    sitename="OpticSimRepeatingStructures.jl",
    format=Documenter.HTML(
        # prettyurls = get(ENV, "CI", nothing) == "true",
        assets=[asset("assets/logo.svg", class=:ico, islocal=true)],
    ),
    modules=[OpticSimRepeatingStructures],
    pages=[
        "Home" => "index.md"
    ]
)

deploydocs(
    repo="github.com/brianguenter/OpticSimRepeatingStructures.jl.git",
    devbranch="main",
    push_preview=true,
)

