using FFTW
using LinearAlgebra
using ProgressBars
using UnPack

function NRCH_2D!(phi0,p,grid,data)
    @unpack N,SO,Tend,dt,A,saveevery = p

    c = Array{ComplexF64}(undef,grid.Nhalf,N,SO)
    c = grid.pfor*phi0
    cdash = similar(c)
    cdashold = similar(c)
    rho = similar(phi0)
    tspan = (2*dt:dt:Tend)
    resize!(data,Int(div(Tend,saveevery))+1)
    data[1] = copy(phi0)
    # push!(data,copy(phi0))
    stride = 2
    println(typeof(c))
    println(typeof(cdash))
    println(typeof(cdashold))
    println(typeof(phi0))
    println(typeof(rho))
    println(typeof(data))
    println(typeof(A))


    NRCH_2D_loop_ETD1!(phi0,rho,c,cdash,cdashold,A,SO,grid)

    for t in ProgressBar(tspan)
        NRCH_2D_loop_ETD2!(phi0,rho,c,cdash,cdashold,A,SO,grid)

        if t%saveevery==0
            data[stride] = copy(phi0)
            stride += 1
            # push!(data,copy(phi0))
        end
    end
end

function NRCH_2D_loop_ETD1!(phi0,rho,c,cdash,cdashold,A,SO,grid)
    computeRho(phi0,rho,A,SO)

    mul!(cdash,grid.pfor,rho)
    cdashold = copy(cdash)

    stepC_ETD1!(c,cdash,grid)

    mul!(phi0,grid.pinv,cdash)
end

function NRCH_2D_loop_ETD2!(phi0,rho,c,cdash,cdashold,A,SO,grid)
    computeRho(phi0,rho,A,SO)

    mul!(cdash,grid.pfor,rho)

    stepC_ETD2!(c,cdash,cdashold,grid)

    mul!(phi0,grid.pinv,cdash)
end

function computeRho(phi,rho,A,SO)
    @inbounds for I in CartesianIndices(rho)
        i,j,k = Tuple(I)
        rho[I] = -phi[I]
        for l=1:SO
            rho[I] += A[k,l]*phi[i,j,l]+phi[i,j,l]^2*phi[I]
        end
    end
end

function stepC_ETD1!(c,cdash,grid)
    @unpack kx2,ky2,kcut2,dt = grid
    cprev = c[1,1,:]
    @inbounds for I in CartesianIndices(c)
        i, j = Tuple(I)
        k2 = kx2[i]+ky2[j]
        if ((kx2[i] < kcut2[1]) | (ky2[j] < kcut2[2]))
            ca = -k2^2
            K1 = exp(ca*dt)
            if abs(ca*dt)<=1.e-5
                K2 = dt+0.5e0*ca*dt^2
            else
                K2 = expm1(ca*dt)/ca
            end
            c[I] = K1 * c[I] - k2 * K2 * cdash[I] 
        else
            c[I] = 0
        end
        cdash[I] = c[I]
    end
    c[1,1,:] = cprev
    cdash[1,1,:] = cprev
end

function stepC_ETD2!(c,cdash,cdashold,grid)
    @unpack kx2,ky2,kcut2,dt = grid
    cprev = c[1,1,:]
    @inbounds for I in CartesianIndices(c)
        i, j = Tuple(I)
        k2 = kx2[i]+ky2[j]
        if ((kx2[i] < kcut2[1]) | (ky2[j] < kcut2[2]))
            ca = -k2^2
            K1 = exp(ca*dt)
            if abs(ca*dt)<=1.e-5
                K2 = 1.5e0*dt+(2.e0/3.e0)*ca*dt^2
                K3 = -0.5e0*dt-(1.e0/6.e0)*ca*dt^2
            else
                K2 = ((1+dt*ca)*expm1(dt*ca)-dt*ca)/(dt*ca^2)
                K3 = (-expm1(dt*ca)+dt*ca)/(dt*ca^2)
            end
            c[I] = K1 * c[I] - k2 * K2 * cdash[I] -k2*K3*cdashold[I]
        else
            c[I] = 0
        end
        cdashold[I] = cdash[I]
        cdash[I] = c[I]
    end
    c[1,1,:] = cprev
    cdash[1,1,:] = cprev
end

function grid2D(p)
    @unpack N,L,SO,dt = p
    Nhalf = div(N,2)+1

    kx= rfftfreq(N,N*(2*pi/L))
    ky = fftfreq(N,N*(2*pi/L))

    kx2 = kx.^2
    ky2 = ky.^2
    kcut = (1/2*maximum(abs.(kx)),1/2*maximum(abs.(ky)))
    kcut2 = kcut.^2
    pfor = plan_rfft(Array{Float64}(undef,N,N,SO),(1,2);flags=FFTW.PATIENT)
    pinv = plan_irfft(Array{ComplexF64}(undef,Nhalf,N,SO),N,(1,2);flags=FFTW.PATIENT)

    return (; N,Nhalf,L,kx2,ky2,kcut2,pfor,pinv,dt)
end

function ASMatrix(a,so)
    # Dont to zeros(Float64,(so,so)) for some reason
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

# Old stuff

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

function stepC_ETDonfly!(c,cdash,grid)
    @unpack kx2,ky2,kcut2,dt = grid
    cprev = c[1,1,:]
    @inbounds for I in CartesianIndices(c)
        i, j = Tuple(I)
        k2 = kx2[i]+ky2[j]
        if ((kx2[i] < kcut2[1]) | (ky2[j] < kcut2[2]))
            ca = -k2^2
            K1 = exp(ca*dt)
            if abs(ca*dt)<=1.e-5
                K2 = dt+0.5e0*ca*dt^2
            else
                K2 = expm1(ca*dt)/ca
            end
            c[I] = K1 * c[I] - k2 * K2 * cdash[I] 
        else
            c[I] = 0
        end
        cdash[I] = c[I]
    end
    c[1,1,:] = cprev
    cdash[1,1,:] = cprev
end


function grid2D_precomp(p)
    @unpack N,L,SO,dt = p

    kx = rfftfreq(N,N*(2*pi/L))
    ky = fftfreq(N,N*(2*pi/L))
    kcut = 2/3*maximum(abs.(kx))
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
            if abs(ca*dt)<=1.e-5
                K2[i,j] = dt+0.5e0*ca*dt^2
            else
                K2[i,j] = expm1(ca*dt)/ca
            end
            K2[i,j] *= -(kx[i]^2+ky[j]^2)
        end
    end
    K2[1,1] = 0
    return (; N,L,kx,ky,ksq,kcut,K1,K2,pfor,pinv)
end

function grid2D_old(p)
    @unpack N,L,SO,dt = p

    kx = rfftfreq(N,N*(2*pi/L))
    ky = fftfreq(N,N*(2*pi/L))
    kcut = 2/3*maximum(abs.(kx))
    pfor = plan_rfft(Array{Float64}(undef,N,N,SO),(1,2);flags=FFTW.PATIENT)
    pinv = plan_irfft(Array{ComplexF64}(undef,div(N,2)+1,N,SO),N,(1,2);flags=FFTW.PATIENT)
    K1 = zeros(div(N,2)+1,N)
    K2 = similar(K1)
    K3 = similar(K1)
    ksq = similar(K1)

    for j in eachindex(ky)
        for i in eachindex(kx)
            ca = -(kx[i]^2+ky[j]^2)^2
            K1[i,j] = exp(ca*dt)
            if abs(ca*dt)<=1.e-5
                K2[i,j] = 1.5e0*dt+(2.e0/3.0e0)*ca*dt^2
                K3[i,j] = -0.5e0*dt-(1.e0/6.e0)*ca*dt^2
            else
                K2[i,j] = ((1+dt*ca)*expm1(dt*ca)-dt*ca)/(dt*ca^2)
                K3[i,j] = (-expm1(dt*ca)+dt*ca)/(dt*ca^2)
            end
            ksq[i,j] = sqrt(kx[i]^2+ky[j]^2)
            K2[i,j] *= -(kx[i]^2+ky[j]^2)
            K3[i,j] *= -(kx[i]^2+ky[j]^2)
        end
    end
    K2[1,1] = 0
    K3[1,1] = 0
    return (; N,L,kx,ky,ksq,kcut,K1,K2,K3,pfor,pinv)
end

