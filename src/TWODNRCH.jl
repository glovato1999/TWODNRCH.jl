module TWODNRCH

include("simulate.jl")
export NRCH_2D!
export grid2D
export ASMatrix

include("render.jl")
export Render

include("readwrite.jl")
export Read
export Write

end
