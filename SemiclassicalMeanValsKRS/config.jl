struct Params

    α        ::Float64                                          # kicking strength -- to match quantum results, this must be negative
    kick     ::Int64                                            # number of kicks
    M0       ::Float64                                          # initial momentum of horizontal lagrangian manifold
    Lₘₐₓ     ::Float64                                          # limit distance between predicted and exact points in final manifold
    Nθ       ::Int64                                            # number of chords used to compute semiclassical mean values
    ħ        ::Vector{Float64}                                  # list of ħ values (each ħ must be chosen as M0/n, n an integer)
    warnings ::Bool                                             # display warnings?
end

struct Detector

    σ   ::Float64     # variance matrix of detector (for now, considered as a multiple of identity -- i.e., circular detector)
    X0  ::SVector{2, Float64}        # center of detector
    R   ::Float64                    # radius of detector to be considered in calculations (for 3 standard deviations, pick σ = 3 × 1/Σ[1,1])
end

# the detector is defined by this "Wigner" function. We assume Σ = something * I
W(X, X0, σ) = inv(2*π*σ^2) * exp( -inv(2σ^2)*norm(X - X0)^2 )

# norm on S¹
d(θ) = min( mod2pi(θ), 2π - mod2pi(θ) )          

# norm on R × S¹, i.e. on the cylinder
D(x) = sqrt( x[1]^2 + d(x[2])^2 )                
