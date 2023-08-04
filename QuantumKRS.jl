""" 
Straightforward implementation of exact quantum evolution to the kicked rotor system. The initial wavefunction in angle representation is

``\\psi_0(\\theta) = \\exp( i \\, m_0 \\theta /\\hbar)/ \\sqrt{2\\pi} \\quad , ``

which amounts to a delta function in action representation. Propagation is performed by iteration:

`` \\vert \\psi_k \\rangle = \\left[ \\exp \\left( i \\, \\alpha \\, \\widehat{m}^2/2\\hbar \\right) \\exp \\left( i \\, \\alpha \\, \\widehat{\\sin \\theta}/\\hbar \\right) \\right] \\vert \\psi_{k-1} \\rangle \\quad , ``

using Fourier grids. 
"""
module QuantumKRS

using StaticArrays: SVector, SMatrix
using FFTW: fft, ifft, plan_ifft, plan_fft, fftshift, fftfreq
using UnPack: @unpack

struct Params

    α      ::Float64                             # kicking strength 
    kick   ::Int64                               # number of kicks
    N      ::Int64                               # size of basis -- note that we always have θ = (0:N-1)*(2π/N), and m the corresponding fftfreqs
    
    M0     ::Float64                             # initial momentum
    ħs     ::Vector{Float64}                     # due to periodicity of wavefunction in angle representation, ħ = M0/n, n being an integrer               
end

struct Detector

    σ  ::Float64                                 # detector's standard deviation
    X0 ::SVector{2, Float64}                     # location of detector
    R  ::Float64                                 # how many detector standard deviations to be included in calculations -- effective radius                             
end

D12(θ1, θ2, ħ, σ, X0) = inv((2π)^1.5 * ħ * σ) * exp( 
                                                    - inv(2σ^2) * ( (θ2+θ1)/2 - X0[2] )^2 +
                                                    - (0.5σ^2/ħ^2) * (θ2-θ1)^2            +
                                                    + (1.0im/ħ) * X0[1] * (θ2-θ1) 
                                                   )
                                                                        
# self-explanatory (A)
function get_grids(N, ħ)

    θ  = (2π/N)    * (0:N-1)
    dθ = abs(θ[2] - θ[1])
    m  = (2π*ħ/dθ) * fftfreq(N)
    dm = abs(m[2] - m[1])

    m, dm, θ, dθ
end

# truncate detector's matrix within θ1[n1:n2] when calculating meal values 
function get_truncation_inds(dθ, X0, R)

    n1 = (Int∘floor)( (X0[2]-R)/dθ )
    n2 = (Int∘floor)( (X0[2]+R)/dθ ) 
    
    n1, n2
end


# wavefunction at final time. Many values of ħ can be passed at once
function wave_at_k(par::Params)

    @unpack α, kick, N, M0, ħs = par
    
    ms   = Vector{Vector{Float64}}([])
    θs   = similar(ms)
    
    ψs_m = Vector{Vector{ComplexF64}}([])
    ψs_θ = similar(ψs_m)
    
    for ħ in ħs
    
        m, dm, θ, dθ = get_grids(N, ħ)

        ψ0  = @. exp( -(1.0im/ħ) * M0 * θ)/sqrt(2π)    # initial state

        ft  = plan_fft(ψ0)
        ift = plan_ifft(ψ0)

        expV = @. exp( -(1.0im/ħ) * α * sin(θ) )
        expK = @. exp( -(1.0im/ħ) * α * m^2/2  )

        for k=1:kick
            ψ0 = ft * ( expK .* ( ift * (expV .* ψ0)))
        end

        ψ_θ = ψ0 
    
        ψ_m = (fftshift∘ifft∘fftshift)(ψ_θ) 
        ψ_m /= sqrt(sum(abs.(ψ_m).^2)*abs(m[2]-m[1]))
        
        push!(ms, fftshift(m))
        push!(θs, θ)
        push!(ψs_m, ψ_m)
        push!(ψs_θ, ψ_θ)
    end
        
    ms, θs, ψs_m, ψs_θ 
end

function mean_vals_at_k(par::Params, mps::Detector)

    @unpack α, kick, N, M0, ħs = par
    @unpack σ, X0, R           = mps
        
    dθ = 2π/N
    
    ms, θs, ψs_m, ψs_θ = wave_at_k(par)
    
    n1, n2 = get_truncation_inds(dθ, X0, R)
    
    Ds = Vector{Float64}([])
    
    for (j,ħ) in enumerate(ħs)

        d12 = [D12(a, b, ħ, σ, X0) for a in θs[j][n1:n2], b in θs[j][n1:n2]]              # truncated detector's matrix
    
        D  = sum( (transpose∘conj)(ψs_θ[j][n1:n2]) * d12 * ψs_θ[j][n1:n2] ) * dθ^2        
        
        imag(D) > 1e-14 && error("Something is wrong: Mean values cannot be complex!")
        
        push!(Ds, real(D))
    end

    Ds
end

end
