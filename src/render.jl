using Plots
using ProgressBars
using StrFormat
using UnPack
using PythonPlot
pythonplot()

function Render(p,data)
    @unpack savename,SO,saveevery = p
    anim = @animate for t in ProgressBar(eachindex(data))
        mini = [minimum(last(data)[j]) for j=1:SO]
        maxi = [maximum(last(data)[j]) for j=1:SO]
        hms = [heatmap(data[t][j][:,:],c=:redsblues,clims=(mini[j],maxi[j]),colorbar_title=f"phi_\%d(j)",aspect_ratio=1) for j=1:SO]
        plot(hms ...,layout=SO)
    end
    mov(anim,"movies/"*savename*".mov",fps=div(60,saveevery),verbose=false,show_msg=false,loop=0)
end