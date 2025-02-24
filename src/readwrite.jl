using HDF5
using ProgressBars
using StrFormat
using UnPack

function Read(datapath)
    data = []
    realpath = "results/"*datapath*".h5"
    param_dic = h5read(realpath,"parameters")
    SO = param_dic["SO"]::Int
    T = param_dic["Tend"]
    N = param_dic["N"]
    saveevery = param_dic["saveevery"]
    tind = [i for i=0:saveevery:T]
    for i in ProgressBar(eachindex(tind))
        phi = []
        for j = 1:SO
            label = f"/phi_\%d(j)/t_\%04d(tind[i])"
            phij = h5read(realpath,label)
            push!(phi,phij)
        end
        push!(data,phi)
    end
    realdata = []
    for i in ProgressBar(eachindex(tind))
        phi = zeros(Float64,N,N,SO)
        for j = 1:SO
            phi[:,:,j] = data[i][j][:,:]
        end
        push!(realdata,phi)
    end
    return param_dic,realdata
end

function Write(p,data)
    @unpack saveevery,Tend,SO,savename = p
    tind = [i for i=0:saveevery:Tend]
    io = h5open("results/"*savename*".h5","w")
    for j = 1:SO
        group_label = f"phi_\%d(j)"
        println("writing ",group_label)
        create_group(io,group_label)
        g = io[group_label]
        for i in ProgressBar(eachindex(data))
            label = f"t_\%04d(tind[i])"
            if haskey(g,label)
                delete_object(g,label)
            end
            g[label]=data[i][:,:,j]
        end
    end

    create_group(io,"parameters")
    g = io["parameters"]

    for key in collect(keys(p))
        strkey = string(key)
        if haskey(g,strkey)
            delete_object(g,strkey)
        end
        g[strkey] = p[key]
    end
    close(io)
end