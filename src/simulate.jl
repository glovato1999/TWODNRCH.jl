using FFTW
using LinearAlgebra
using ProgressBars
using UnPack

function NRCH_2D!(phi0,p,grid,data)
    @unpack N,SO,Tend,dt,A,saveevery = p

    Nhalf = div(N,2)+1
    c = Array{ComplexF64}(undef,Nhalf,N,SO)
    c = grid.pfor*phi0
    cdash = similar(c)
    rho = similar(phi0)
    tspan = [i for i=0:dt:Tend]
    resize!(data,Int(div(Tend,saveevery))+1)
    stride = 1
    for t in ProgressBar(tspan)
        NRCH_2D_loop!(phi0,rho,c,cdash,A,SO,grid)

        if t%saveevery==0
            data[stride] = copy(phi0)
            stride += 1
        end
    end
end

function NRCH_2D_loop!(phi0,rho,c,cdash,A,SO,grid)
    computeRho(phi0,rho,A,SO)

    mul!(cdash,grid.pfor,rho)

    stepC!(c,cdash,grid)
    # stepC_ETDonfly!(c,cdash,grid,dt)

    mul!(phi0,grid.pinv,cdash)
end

function computeRho(phi,rho,A,SO)
    @inbounds for I in CartesianIndices(rho)
        i,j,k = Tuple(I)
        rho[I] = -phi[I]+phi[I]^3
        for l=1:SO
            rho[I] += A[k,l]*phi[i,j,l]
        end
    end
end

function stepC!(c,cdash,grid)
    @unpack K1,K2,ksq,kcut = grid
    @inbounds for I in CartesianIndices(c)
        i, j = Tuple(I)
        if ksq[i,j] < kcut
            c[I] = K1[i,j] * c[I] + K2[i,j] * cdash[I] 
        else
            c[I] = 0
        end
        cdash[I] = c[I]
    end
end

function stepC_ETDonfly!(c,cdash,grid,dt)
    @unpack kx2,ky2,kcut2 = grid
    cprev = c[1,1,:]
    @inbounds for I in CartesianIndices(c)
        i, j = Tuple(I)
        k2 = kx2[i]+ky2[j]
        if k2 < kcut2
            ca = -k2^2
            K1 = exp(ca*dt)
            K2 = expm1(ca*dt)/k2
            c[I] = K1 * c[I] + K2 * cdash[I] 
        else
            c[I] = 0
        end
        cdash[I] = c[I]
    end
    c[1,1,:] = cprev
    cdash[1,1,:] = cprev
end

function grid2D(p)
    @unpack N,L,SO,dt = p

    kx = rfftfreq(N,N*(2*pi/L))
    ky = fftfreq(N,N*(2*pi/L))
    kcut = 1/2*maximum(abs.(kx))
    pfor = plan_rfft(Array{Float64}(undef,N,N,SO),(1,2);flags=FFTW.PATIENT)
    pinv = plan_irfft(Array{ComplexF64}(undef,div(N,2)+1,N,SO),N,(1,2);flags=FFTW.PATIENT)
    K1 = zeros(div(N,2)+1,N)
    K2 = similar(K1)
    ksq = similar(K1)

    for j in eachindex(ky)
        for i in eachindex(kx)
            ca = -(kx[i]^2+ky[j]^2)^2
            ksq[i,j] = sqrt(kx[i]^2+ky[j]^2)
            K1[i,j] = exp(ca*dt)
            K2[i,j] = -ksq[i,j]*(exp(ca*dt)-1)/ca
        end
    end
    K2[1,1] = 0
    return (; N,L,kx,ky,ksq,kcut,K1,K2,pfor,pinv)
end

function grid2D_ETDonfly(p)
    @unpack N,L,SO = p

    kx2 = rfftfreq(N,N*(2*pi/L)).^2
    ky2 = fftfreq(N,N*(2*pi/L)).^2
    kcut2 = (1/2)^2*maximum(abs.(kx2))
    pfor = plan_rfft(Array{Float64}(undef,N,N,SO),(1,2);flags=FFTW.PATIENT)
    pinv = plan_irfft(Array{ComplexF64}(undef,div(N,2)+1,N,SO),N,(1,2);flags=FFTW.PATIENT)

    return (; N,L,kx2,ky2,kcut2,pfor,pinv)
end

