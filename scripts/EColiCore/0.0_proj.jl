# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
using RegulationImageMC_2024
using ProjFlows
using MetX
using MetX.MetXGEMs
using MetX.MetXEP
using Bloberias

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
#  Project
# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.

PROJ = Project0(RegulationImageMC_2024)
globproj!(PROJ)

# TODO: implement version into project0
# TODO: implement a `versioninfo`-like method for Projects to summarize 
PROJVER = v"0.2.0-ecoli.core"

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
include("0.99_proj.base.jl")

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# Bloberia
B = Bloberia(_blobsdir())
G = blob!(B, "sim.globals")
mkpath(B)

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
return nothing