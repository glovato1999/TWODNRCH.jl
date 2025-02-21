module TWODNRCH

include("simulate.jl")
export NRCH_2D!
export grid2D
export grid2D_ETDonfly
export ASMatrix

include("render.jl")
export Render
export Render_Amp

include("readwrite.jl")
export Read
export Write

end
