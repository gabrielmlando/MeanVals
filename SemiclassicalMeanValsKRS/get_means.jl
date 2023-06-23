# calculates SMV means using ħ as a parameter 
function get_means_smv(par::Params, dec::Detector)
    
    @unpack α, kick, M0, Lₘₐₓ, Nθ, ħ, warnings = par
    @unpack R, Σ, X0                           = dec
    
    θ0s, Lk                    = get_Lk(α, kick, M0, Lₘₐₓ, X0, R, warnings)                  # obtain sections of Lk falling inside the detector
    θ_fins, fins, θ_fils, fils = filament_or_finger(θ0s, Lk)                                 # check which are filaments and which are fingers                    
    
    sort_from_left_to_right!(M0, θ_fils, fils)                                               # sort filaments from left to right, then up to down
                                  
    θ_fils, fils               = sort_from_up_to_down(θ_fils, fils)  
    
    classical_result           = get_means_cla(α, kick, M0, X0, Σ, Nθ, θ0s, Lk)              # get constant term <O>_cl
    semiclassical_result       = get_means_osc(α, kick, M0, X0, Σ, Nθ, θ_fils, fils, ħ)      # get oscillatory term <O>_sc

    if warnings
        println("Number of filaments inside the detector: ", length(θ_fils))
        println("Number of fingers inside the detector  : ", length(θ_fins))
    end
    
    classical_result .+ semiclassical_result                                                 # return <O>_cl + <O>_sc
end

