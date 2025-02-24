using Plots
using ProgressBars
using StrFormat
using UnPack

function Render(p,data)
    @unpack savename,SO,saveevery = p
    anim = @animate for t in ProgressBar(eachindex(data))
        mini = [minimum(last(data)[:,:,j]) for j=1:SO]
        maxi = [maximum(last(data)[:,:,j]) for j=1:SO]
        dm = maxi .- mini
        hms = [heatmap(data[t][:,:,j],c=:redsblues,clims=(mini[j],maxi[j]),title=f"phi_\%d(j)",aspect_ratio=1,axis=([],false),cbar=false) for j=1:SO]
        plot(hms ...,layout=SO)
    end
    mov(anim,"movies/"*savename*".mov",fps=div(60,saveevery),verbose=false,show_msg=false,loop=0)
end

function Render_Amp(p,data)
    @unpack savename,SO,saveevery,N = p
    renderdata = []
    for i in ProgressBar(eachindex(data))
        phimod = zeros(Float64,N,N)
        for J in CartesianIndices(data[i])
            j,k,l = Tuple(J)
            phimod[j,k] += data[i][J].^2
        end
        push!(renderdata,copy(phimod))
    end
    mini = minimum(last(renderdata))
    maxi = maximum(last(renderdata))
    anim = @animate for i in ProgressBar(eachindex(renderdata))
        plot(heatmap(renderdata[i],c=:greys,clims=(mini,maxi)))
    end
    mov(anim,"movies/"*savename*".mov",fps=div(60,saveevery),verbose=false,show_msg=false,loop=0)
end