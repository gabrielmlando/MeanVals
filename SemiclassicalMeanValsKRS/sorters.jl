function filament_or_finger(θ0s, Lk)
    
    θ_fins = Vector{Vector{Float64}}([]); fins = Vector{Vector{SVector}}([])
    θ_fils = Vector{Vector{Float64}}([]); fils = Vector{Vector{SVector}}([])

    for k=1:length(Lk)

        if sign( last(Lk[k][2] - Lk[k][1])) != sign( last(Lk[k][end] - Lk[k][end-1]))
            push!(θ_fins, θ0s[k])
            push!(fins, Lk[k])
        else
            push!(θ_fils, θ0s[k])
            push!(fils, Lk[k])
        end
    end
    
    θ_fins, fins, θ_fils, fils
end

function sort_from_left_to_right!(M0, θ_fils, fils)
    
    for k=1:length(fils)
        if D(fils[k][1]-[M0,0.0]) < D(fils[k][end]-[M0,0.0]) 
            fils[k]   = reverse(fils[k])
            θ_fils[k] = reverse(θ_fils[k])
        end
    end
end

function sort_from_up_to_down(θ_fils, fils)                       # input here needs to have been sorted by sort_from_left_to_right! 
    
    srt      = sortperm(norm.(first.(fils) - last.(fils)))        # creates table of indices ordered according to norm of filaments, so they will be 
    fils_srt = fils[srt]                                          # interwined between the two hemispheres (smallest to largest). Largest index is [end].
    θs_srt   = θ_fils[srt]
    
    upper_fils = Vector{AbstractVector}([])
    lower_fils = Vector{AbstractVector}([])
    
    upper_θs   = Vector{AbstractVector}([])
    lower_θs   = Vector{AbstractVector}([])
      
    for k=1:length(fils_srt)
    
        m_diff = first(first(fils_srt[k]) - first(fils_srt[end])) # this is the differente in m's between filament [k] and the largest filament, [end]
                                                                  # Note that, by obvious reasons, I am assuming the largest filament lies near the equator
        if m_diff < 0.0                                           # If m[k] - m[end] is smaller than 0, it means filament [k] is below the equator
            
            push!(lower_fils, fils_srt[k]) 
            push!(lower_θs,   θs_srt[k])
            
        elseif m_diff >= 0.0                                      # If m[k] - m[end] is larger than 0, it means filament [k] is above the equator or k=end
            
            push!(upper_fils, fils_srt[k])
            push!(upper_θs,   θs_srt[k])
        end
    end
    
    vcat(upper_θs, reverse(lower_θs)), vcat(upper_fils, reverse(lower_fils))
end

