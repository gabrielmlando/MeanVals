""" 
Straightforward implementation of the Herman-Kluk wavefunction to the kicked rotor system in module `QuantumKRS`. Classical propagation follows

`` m_k = m_{k-1} - \\alpha \\, \\cos \\theta_{k-1} \\, \\quad \\theta_k = \\theta_{k-1} + \\alpha \\, m_k \\quad (\\rm{mod} \\, 2\\pi) \\quad .``

For latest update, see https://github.com/gabrielmlando/MeanVals
"""
module HermanKlukKRS

using StaticArrays: SMatrix
using UnPack: @unpack
using Distributions: Normal, Uniform
using LinearAlgebra: I, mul!
using FFTW: fft, fftfreq, fftshift

struct Params

    α      ::Float64              # kicking strength         
    kick   ::Int64                # number of kicks
    N      ::Int64                # size of basis -- note that we always have θ = (0:N-1)*(2π/N), and m the corresponding fftfreqs (see (A) below)
    
    M0     ::Float64              # initial momentum
    ħ      ::Float64              # due to periodicity of wavefunction in angle representation, ħ = M0/n, n being an integrer
    
    γ      ::Float64              # coherent state width
    ntrajs ::Int64                # number of trajectories for Monte Carlo
end

mutable struct Variables
    
    mₖ₋₁::Float64                 # final action
    θₖ₋₁::Float64                 # final angle
    sₖ₋₁::Float64                 # action
    Aₖ₋₁::SMatrix{2, 2, Float64}  # monodromy matrix
    rₖ₋₁::ComplexF64              # HK prefactor
    μₖ₋₁::Int64                   # exp(1.0im*π*μ) = (-1)^μ = ±1, where μ = number of pre-factor branch crossings 
    
    Variables(m₀, θ₀) = new(m₀, θ₀, 0.0, SMatrix{2,2}(1.0I), 1.0+0.0im, 1)
end

# coherent state on the circle in action representation
g(y, m, θ, γ, ħ) = (π*γ*ħ)^(-1/4) * exp( -inv(2.0γ*ħ)*(y-m)^2 - (1.0im/ħ)*y*θ )

# unnormalized MCMC estimator -- the posterior is the real part of coherent state above, i.e. exp( -inv(2.0σ)*(y-m)^2 ), σ = γ*ħ
w(m, θ, γ, ħ) = exp( (1.0im/ħ) * m * θ)

# self-explanatory (A)
function get_grids(N, ħ)

    θ  = (2π/N)    * (0:N-1)
    dθ = abs(θ[2] - θ[1])
    m  = (2π*ħ/dθ) * fftshift(fftfreq(N))
    dm = abs(m[2] - m[1])

    m, dm, θ, dθ
end

# one iteration of KRS
function step!(par::Params, var::Variables)  
    
    @unpack α, γ, ħ = par
    
    @unpack mₖ₋₁, θₖ₋₁, sₖ₋₁, Aₖ₋₁,  rₖ₋₁, μₖ₋₁ = var
    
    # map
    mₖ = mₖ₋₁ - α*cos(θₖ₋₁)        
    θₖ = mod2pi(θₖ₋₁ + α*mₖ)      
    
    # action
    sₖ = sₖ₋₁ + α*(0.5*mₖ^2 - sin(θₖ₋₁))
    
    # monodromy: [∂mₖ/∂m₀  ∂mₖ/∂θ₀ ; ∂θₖ/∂m₀  ∂θₖ/∂θ₀]
    Atemp = SMatrix{2,2}([   1.0          α*sin(θₖ₋₁)      ;      
                              α       1.0 + α^2*sin(θₖ₋₁)  ])      
    Aₖ = mul!(similar(Atemp), Atemp, Aₖ₋₁)
    
    # pre-factor: (1/2)( ∂mₖ/∂m₀ + ∂θₖ/∂θ₀ ) + (i/2)( (1/γ)*∂mₖ/∂θ₀ - γ*∂θₖ/∂m₀ )
    rₖ = 0.5*(Aₖ[1,1] + Aₖ[2,2]) + 0.5im*( inv(γ)*Aₖ[1,2] - (γ)*Aₖ[2,1] ) 
    
    # branch changes
    μₖ = (real(rₖ) < 0.0) && (imag(rₖ)*imag(rₖ₋₁) < 0.0) ? -μₖ₋₁ : μₖ₋₁
   
    # update variables
    var.mₖ₋₁ = mₖ
    var.θₖ₋₁ = θₖ
    var.sₖ₋₁ = sₖ
    var.Aₖ₋₁ = Aₖ
    var.rₖ₋₁ = rₖ
    var.μₖ₋₁ = μₖ
end

# iterate step! to produce orbit + action + monodromy + pre-factor + branch changes (monodromy is not returned)
function flow!(par::Params, m₀, θ₀)
    
    var = Variables(m₀, θ₀)
    
    for _ = 1:par.kick 
        step!(par, var)
    end
    
    var.mₖ₋₁, var.θₖ₋₁, var.sₖ₋₁, var.rₖ₋₁, var.μₖ₋₁
end
  
# wavefunction at kick=k
function wave_at_k(par::Params)

    @unpack α, kick, N, γ, M0, ħ, ntrajs = par
    
    m, dm, θ, dθ = get_grids(N, ħ)
     
    mdist = Normal(M0, sqrt(γ*ħ)) 
    θdist = Uniform(0, 2π)
    
    # wavefunction computed in action representation
    ψ_m = zeros(ComplexF64, N) 

    for _=1:ntrajs 

        m₀ = rand(mdist)
        θ₀ = rand(θdist)
        
        mₖ, θₖ, sₖ, rₖ, μₖ = flow!(par, m₀, θ₀)
        
        for j=1:N
            
            ψ_m[j] += μₖ * sqrt(rₖ) *                       # pre-factor  
                      w(M0, θ₀, γ, ħ) *                     # MCMC estimator
                      g(m[j], mₖ, θₖ, γ, ħ) *               # dynamic coherent state in angle representation
                      exp(1.0im*sₖ/ħ)                       # exponential of classical action      
        end    
    end
    
    # renormalize
    ψ_m /= sqrt(sum(abs.(ψ_m).^2)*dm) 
    
    # wavefunction in angle representation + renormalization
    ψ_θ =  (fftshift∘fft∘fftshift)(ψ_m)
    ψ_θ /= sqrt(sum(abs.(ψ_θ).^2)*dθ) 
    
    m, θ, ψ_m, ψ_θ 
end 

end