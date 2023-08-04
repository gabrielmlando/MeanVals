# computes oscillatory part of SMV mean values
function get_means_osc(α, kick, m0, X0, σ, Nθ, θ_fils, fils, ħs)
    
    N = length(fils)

    result = zeros(Float64, length(ħs))
    
    Threads.@threads for i=1:N-1
        Threads.@threads for j=i+1:N
                  
            θ1 = θ_fils[i][1]
            θ2 = θ_fils[j][1]
            
            # parametrization increment is constant for each filament and doesn't need to be updated inside k loop
            dθ1 = (θ_fils[i][end]-θ_fils[i][1]) / (Nθ - 1)
            dθ2 = (θ_fils[j][end]-θ_fils[j][1]) / (Nθ - 1)
            
            # start integral 
            ∫dη = zeros(Float64, length(ħs))
            
            for k=1:Nθ
                
                m1′, θ1′, s1, A1, μ1 = flow_ρ!(α, kick, SVector{2}(m0, θ1))
                m2′, θ2′, s2, A2, μ2 = flow_ρ!(α, kick, SVector{2}(m0, θ2))

                ###################################################
                ##################### AREAS #######################
                ###################################################

                s▱ = 0.5*(m2′ + m1′)*(θ2′ - θ1′)                         
                s□ = m0*(θ2 - θ1)                                      
      
                Δs = (s2 - s1) - s▱ + s□
                 
                # center between filaments
                η  = 0.5*([m1′, θ1′] + [m2′, θ2′]) 
   
                for (r,ħ) in enumerate(ħs)
                    
                    ∫dη[r] += W(η, X0, σ) * cos(Δs/ħ + 0.5π*(μ2 - μ1) ) * (sqrt∘abs)(dθ1*dθ2)
                end
                
                ###################################################
                #################### UPDATE #######################
                ###################################################
                
                θ1   += dθ1
                θ2   += dθ2
            end
            
            for r=1:length(ħs)
                result[r] += ∫dη[r]/π
            end
        end
    end
    
    result
end
