using TWODNRCH
using HDF5
using Test
using ProgressBars
using StrFormat
using Plots
using Statistics
using LinearAlgebra
using FFTW

@testset "TWODNRCH.jl" begin
    gr()
    N = 256
    SO = 2
    Tend = 1e2
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

    path = "test/real_NRCH_ETD1.h5"
    refdata_1 = []
    tind = (0:saveevery:Tend)
    for i in ProgressBar(tind)
        phi = zeros(Float64,N,N,SO)
        phi[:,:,1] = h5read(path,f"/phi_1/t_\%04d(i)")
        phi[:,:,2] = h5read(path,f"/phi_2/t_\%04d(i)")
        push!(refdata_1,copy(phi))
    end

    path = "test/2D_NRCH.h5"
    refdata_2 = []
    tind = (0:saveevery:Tend)
    for i in ProgressBar(tind)
        phi = zeros(Float64,N,N,SO)
        phi[:,:,1] = h5read(path,f"/phi_1/t_\%04d(i)")
        phi[:,:,2] = h5read(path,f"/phi_2/t_\%04d(i)")
        push!(refdata_2,copy(phi))
    end 

    phi0 = copy(refdata_1[1])

    tpoints_1,data_1 = NRCH_2D(phi0,p,1)

    phi0 = copy(refdata_2[1])
    tpoints_2,data_2 = NRCH_2D(phi0,p,2)

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

    end

    p = plot(tind[1:end],diffmeasure_ref1my1[1:end],label="ETD1")
    plot!(p,tind[1:end],diffmeasure_ref2my2[1:end],label="ETD2")
    savefig(p,"test/data_error.svg")
end


function tw(t,r,q,alpha)
    theta = (-alpha*dot(q,q)*t-dot(q,r))
    rho = sqrt(1-dot(q,q))
    return rho*[cos(theta); sin(theta)]
end

@testset "TravelingWave" begin
    gr()
    N = L = 32
    SO = 2
    alpha = 0.2
    A = alpha*[0. -1.; 1. 0]
    p = (;
        N,
        L,
        SO,
        dt = 5.e-2,
        Tend = 1.e2,
        saveevery = 4,
        A
    )
    
    dx = L/(N-1)
    x=y=[i*dx for i=0:N-1]
    kx=ky = rfftfreq(N,N*2*pi/L)
    q = [kx[3] 0]
    
    phi0 = zeros(Float64,N,N,SO)
    
    for i=1:N
        for j=1:N
            r = [x[i] y[j]]
            phi0[i,j,:] .= tw(0,r,q,alpha)
        end
    end
    
    tpoints,data_1 = NRCH_2D(phi0,p,1)
    tpoints,data_2 = NRCH_2D(phi0,p,2)
    
    
    data_tw = Array{Float64,3}[]
    for ti in eachindex(tpoints)
        phisol = zeros(Float64,N,N,SO)
        for i=1:N
            for j=1:N
                r = [x[i] y[j]]
                phisol[i,j,:] = tw(tpoints[ti],r,q,alpha)
            end
        end
        push!(data_tw,copy(phisol))
    end
    
    onepoint_1 = []
    onepoint_2 = []
    onepoint_ana = []
    for ti in eachindex(data_1)
        push!(onepoint_1,copy(data_1[ti][1,1,1]))
        push!(onepoint_2,copy(data_2[ti][1,1,1]))
        push!(onepoint_ana,copy(data_tw[ti][1,1,1]))
    end
    pl2 = plot(tpoints,onepoint_1,label="ETD_1")
    plot!(pl2,tpoints,onepoint_2,label="ETD_2")
    plot!(pl2,tpoints,onepoint_ana,label="analytical")
    savefig(pl2,"test/onepoint.svg")
end

