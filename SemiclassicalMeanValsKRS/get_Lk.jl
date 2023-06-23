#= Computes the evolved Lagrangian manifold Lk, but stores only the portions falling inside the detector. =#
function get_Lk(α, kick, M0, Lₘₐₓ, X0, R, warnings)

    θ0s_temp  = Vector{Float64}([])         # will be updated as many times as there are filaments      
    θ0s       = Vector{Vector{Float64}}([])  # an array of arrays, each representing the parametrization of a filament
 
    Lk_temp = Vector{SVector}([])
    Lk      = Vector{Vector{SVector}}([])
    
    L  = Lₘₐₓ                              # maximum distance between predictions using linearized and true dynamics 

    θ0 = 0.0; θ0_loop = copy(θ0)
    xk_old, Ak_old = flow_μ!(α, kick, SVector{2}(M0, θ0)) 
    
    D(xk_old - X0) < R && error("Initial point fell inside the detector.")

    m_tot = 0; m_rej = 0                   # m counts failures
    while θ0 < 2π 
        
        v  = Ak_old * [0.0, 1.0]           # v ∈ T(S¹ × R)
        
        δθ = L/norm(v)                     # Tangent space is R², so R² norm is used instead of cylinder norm
        
        xk_pred = xk_old + v * δθ          # xk predicted by linearization
        
        xk_new, Ak_new = flow_μ!(α, kick, SVector{2}(M0, θ0+δθ)) 
        
        γ = D(xk_pred - xk_new)            # distance between predicted xk and computed xk

        if γ < 0.25Lₘₐₓ                    # if smaller than 25% of Lₘₐₓ, point is accepted
            
            L      = Lₘₐₓ                  # keep Lₘₐₓ as is
            θ0    += δθ                    # move forward in the parametrization angle
            
            if D(xk_new - X0) <= R         # if point falls inside detector
                push!(Lk_temp, xk_new)     # push point to temporary filament
                push!(θ0s_temp , θ0)       # together with parametrization angle
            end
            
            if D(xk_new - X0) > R && D(xk_old - X0) <= R   # if the point falls outside of the detector but comes from the inside,
                push!(Lk, Lk_temp)                         # the filament is considered to be finished, so push temporary filament to Lk
                push!(θ0s , θ0s_temp)                      # together with parametrization angle
                
                
                length(Lk_temp) > 3 || begin pop!(Lk); pop!(θ0s); # if fillament is too short, code will complain that this might trigger errors
                                             if warnings
                                                @warn ("Some portions of Lk are not covered. Smaller Lₘₐₓ should fix this.") 
                                             end
                                       end
                
                Lk_temp = Vector{SVector}([])                 # reinitiate filaments
                θ0s_temp = Vector{Float64}([])                # and angles
            end
            
            xk_old = xk_new; Ak_old = Ak_new # also, since the point is accepted, we must update the point and the monodromy matrix

        else           # if the point is declined                              
            L /= 2     # try again with half the initial step
            m_rej += 1 # counts the number of rejected points
        end
        
        m_tot += 1
    end

    if warnings
        println("Percentage of rejected trajectories: ", 100*round(m_rej/m_tot, digits=4), "%")
    end
    
    θ0s, Lk
end;
