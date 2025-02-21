using TWODNRCH
using HDF5
using Test
using ProgressBars
using StrFormat
using Plots
using Statistics

@testset "TWODNRCH.jl" begin
    gr()
    N = 256
    SO = 2
    Tend = 1e2
    # Tend = 4
    dt = 5.e-2
    savename = ""
    a = ones(Float64,div(SO*(SO-1),2))
    alpha = -0.2
    NonR = ASMatrix(a,SO)
    saveevery = 4
    p = (;
        N,
        L=N,
        SO,
        Tend,
        dt,
        savename,
        A=alpha*NonR,
        saveevery
    )
    g = grid2D(p)
    println(typeof(g))

    path = "test/2D_NRCH_ETD1.h5"
    refdata_1 = []
    tind = (0:saveevery:Tend)
    for i in ProgressBar(tind)
        phi = zeros(Float64,N,N,SO)
        phi[:,:,1] = h5read(path,f"/phi_1/t_\%04d(i)")
        phi[:,:,2] = h5read(path,f"/phi_2/t_\%04d(i)")
        push!(refdata_1,copy(phi))
    end

    path = "test/2D_NRCH_ETD1.h5"
    refdata_2 = []
    tind = (0:saveevery:Tend)
    for i in ProgressBar(tind)
        phi = zeros(Float64,N,N,SO)
        phi[:,:,1] = h5read(path,f"/phi_1/t_\%04d(i)")
        phi[:,:,2] = h5read(path,f"/phi_2/t_\%04d(i)")
        push!(refdata_2,copy(phi))
    end 

    phi0 = zeros(Float64,N,N,SO)
    phi0 = copy(refdata_1[1])

    data_1 = []
    NRCH_2D!(phi0,p,g,data_1)
    println(typeof(data_1))

    println(typeof(g))
    phi0 = copy(refdata_2[1])
    data_2 = []
    NRCH_2D!(phi0,p,g,data_2)
    println(typeof(data_2))

    diffmeasure_ref12 = []
    diffmeasure_ref1my1 = []
    diffmeasure_ref2my2 = []
    diffmeasure_my12 = []
    for i in eachindex(data_2)
        amp_ref1 = refdata_1[i][:,:,1]
        amp_ref2 = refdata_2[i][:,:,1]
        amp_my1 = data_1[i][:,:,1]
        amp_my2 = data_2[i][:,:,1]
        diff_ref12 = abs.(amp_ref1-amp_ref2)
        diff_ref1my1 = abs.(amp_ref1-amp_my1)
        diff_ref2my2 = abs.(amp_ref2-amp_my2)
        diff_my12 = abs.(amp_my1-amp_my2)
        push!(diffmeasure_my12,mean(diff_my12))
        push!(diffmeasure_ref12,mean(diff_ref12))
        push!(diffmeasure_ref1my1,mean(diff_ref1my1))
        push!(diffmeasure_ref2my2,mean(diff_ref2my2))

        # println("comp to ETD1 ",mean(diff_ETD1),", comp to ETD2 ",mean(diff_ETD2))
        # pms = (heatmap(diff_ETD1,clims=(0,0.01)),heatmap(diff_ETD2,clims=(0,0.01)))

        # display(plot(pms...,layout=2))
    end
    # println(diffmeasure_my12[1])
    # p = plot(tind[2:end],diffmeasure_my12[2:end],label="my12",xscale=:log10,yscale=:log10)
    p = plot(tind[2:end],diffmeasure_my12[2:end],label="my12")
    # plot!(p,tind,diffmeasure_ref12,label="ref12")
    # plot!(p,tind[2:end],diffmeasure_ref1my1[2:end],label="ref1my1")
    # plot!(p,tind[2:end],diffmeasure_ref2my2[2:end],label="ref2my2")
    display(p)
end
