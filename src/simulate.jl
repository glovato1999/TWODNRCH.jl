using FFTW
using LinearAlgebra
using ProgressBars
using UnPack
using StrFormat

function grid2D(p)
    @unpack N,L,SO = p

    Nhalf = div(N,2)+1

    kx = rfftfreq(N,N*2*pi/L)
    ky = fftfreq(N,N*2*pi/L)
    kx2 = kx.^2
    ky2 = ky.^2
    kcut = 2/3*maximum(abs.(kx))
    kcut2 = kcut^2

    pfor = plan_rfft(Array{Float64}(undef,N,N,SO),(1,2),flags=FFTW.PATIENT)
    pinv = plan_irfft(Array{ComplexF64}(undef,Nhalf,N,SO),N,(1,2),flags=FFTW.PATIENT)

    return (; kx2,ky2,kcut2,pfor,pinv,Nhalf)
end

function NRCH_2D(phi0,parameters,ETD)
    g = grid2D(parameters)
    @unpack N,SO,dt,Tend,saveevery,A = parameters
    @unpack Nhalf,kx2,ky2,kcut2,pfor,pinv = g

    phik = zeros(ComplexF64,Nhalf,N,SO)
    nonlin = zeros(Float64,N,N,SO)
    nonlink = zeros(ComplexF64,Nhalf,N,SO)
    phi = copy(phi0)

    t = (dt:dt:Tend)

    data = Array{Float64,3}[]
    tpoints = Float64[]
    resize!(data,Int(div(Tend,saveevery))+1)
    resize!(tpoints,Int(div(Tend,saveevery))+1)
    data[1] = copy(phi)
    tpoints[1] = 0
    stride = 2

    mul!(phik,pfor,phi)

    if (ETD == 1)
        println("ETD_1")
        @time begin
            for ti in eachindex(t)
                perc = ti/size(t)[1]*100
                print(f"\%.1f(perc)% \%.1f(t[ti])/\%.1f(Tend)s\r")

                computeNonlin!(nonlin,phi,A,SO)

                mul!(nonlink,pfor,nonlin)

                step_ETD1!(phik,nonlink,kx2,ky2,kcut2,dt)

                mul!(phi,pinv,nonlink)

                if t[ti]%saveevery == 0
                    data[stride] = copy(phi)
                    tpoints[stride] = t[ti]
                    stride += 1

                end
            end
        end
    elseif (ETD == 2)
        println("ETD_2")
        t = (2*dt:dt:Tend)
        nonlink_old = zeros(ComplexF64,Nhalf,N,SO)
        @time begin
            computeNonlin!(nonlin,phi,A,SO)
            mul!(nonlink,pfor,nonlin)
            copyto!(nonlink_old,nonlink)
            step_ETD1!(phik,nonlink,kx2,ky2,kcut2,dt)
            mul!(phi,pinv,nonlink)

            for ti in eachindex(t)
                perc = ti/size(t)[1]*100
                print(f"\%.1f(perc)% \%.1f(t[ti])/\%.1f(Tend)s\r")

                computeNonlin!(nonlin,phi,A,SO)
                mul!(nonlink,pfor,nonlin)

                step_ETD2!(phik,nonlink,nonlink_old,kx2,ky2,kcut2,dt)

                mul!(phi,pinv,nonlink)

                if t[ti]%saveevery == 0
                    data[stride] = copy(phi)
                    tpoints[stride] = t[ti]
                    stride += 1
                end
            end
        end
    end
    return tpoints,data
end

function computeNonlin!(nonlin,phi,A,SO)
    for I in CartesianIndices(nonlin)
        i,j,k = Tuple(I)
        nonlin[I] = -phi[I]
        for l=1:SO
            nonlin[I] += phi[i,j,l]^2*phi[I]
            nonlin[I] += A[k,l]*phi[i,j,l]
        end
    end
    nothing
end

function step_ETD1!(phik,nonlink,kx2,ky2,kcut2,dt)
    for I in CartesianIndices(phik)
        i,j,k = Tuple(I)
        if (kx2[i] < kcut2) & (ky2[j] < kcut2)
            k2 = kx2[i]+ky2[j]
            ca = -k2^2
            K1 = exp(ca*dt)
            if (abs(ca*dt) <= 1.e-5)
                K2 = dt + 0.5e0 * ca * dt^2
            else
                K2 = expm1(ca*dt)/ca
            end
            phik[I] = K1*phik[I] - k2*K2*nonlink[I]
        else
            phik[I] = 0
        end
        nonlink[I] = phik[I]
    end
    nothing
end

function step_ETD2!(phik,nonlink,nonlink_old,kx2,ky2,kcut2,dt)
    for I in CartesianIndices(phik)
        i,j,k = Tuple(I)
        if (kx2[i] < kcut2) & (ky2[j] < kcut2)
            k2 = kx2[i]+ky2[j]
            ca = -k2^2 
            K1 = exp(dt*ca)
            if (abs(dt*ca) <= 1.e-5)
                K2 = 1.5e0 * dt + (2.e0 / 3.e0) * ca * dt^2
                K3 = -0.5e0 * dt - (1.e0 / 6.e0) * ca * dt^2
            else
                K2 = ((1 + dt * ca) * expm1(dt * ca) - dt * ca) / (dt*ca^2)
                K3 = (-expm1(dt * ca) + dt * ca) / (dt * ca^2)
            end
            phik[I] = K1 * phik[I] - k2 * K2 * nonlink[I] - k2 * K3 * nonlink_old[I]
        else
            phik[I] = 0
        end
        nonlink_old[I] = nonlink[I]
        nonlink[I] = phik[I]
    end
end

function ASMatrix(a,so)
    A = zeros(Float64,so,so)
    signum = -1
    strife = 1
    for i=1:so
        for j=i+1:so
            A[i,j]=signum*a[strife]
            A[j,i]=-signum*a[strife]
            signum*=-1
            strife+=1
        end
    end
    return A
end
