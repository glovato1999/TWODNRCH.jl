module TWODNRCH


include("simulate.jl")
export NRCH_2D!
export NRCH_2D
export grid2D
export grid2D_ETDonfly
export ASMatrix
export tw

include("render.jl")
export Render
export Render_Amp

include("readwrite.jl")
export Read
export Write

end
