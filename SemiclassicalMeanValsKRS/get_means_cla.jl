# calculate constant (classical) contribution to smv means using the truncated Wigner approximation
function get_means_cla(α, kick, M0, X0, Σ, Nθ, θ0s, Lk)

    N = length(θ0s)
    
    result = 0.0
    for i=1:N

        dθ = abs(θ0s[i][end]-θ0s[i][1]) / (Nθ - 1)
        θₖ = θ0s[i][1]

        # start integral 
        ∫dη = 0.0 
        for k=1:Nθ
            
            xₖ′ = flow!(α, kick, SVector{2}(M0, θₖ))
        
            ∫dη += W(xₖ′, X0, Σ)*dθ 
        
            θₖ += dθ
        end
        result += ∫dη/2π
    end
    result
end
