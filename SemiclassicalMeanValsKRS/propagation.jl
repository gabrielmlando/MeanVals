function ∇H!(α, x)   
    # map
    M = x[1] - α*cos(x[2])        
    Θ = mod2pi(x[2] + α*M)      
    
    SVector{2}(M, Θ)
end

function flow!(α, kick, x)

    for _=1:kick  
        x = ∇H!(α, x)
    end
    x
end


# dynamics returning final point and monodromy. Used in get_Lk
function ∇H_μ!(α, x, A_old)   
    # map
    M = x[1] - α*cos(x[2])        
    Θ = mod2pi(x[2] + α*M)      
    
    # monodromy
    A_new = SMatrix{2,2}([   1.0        α*sin(x[2])      ;    
                              α   1.0 + α^2*sin(x[2])]    )
    
    A_new = mul!(similar(A_old), A_new, A_old)
    
    SVector{2}(M, Θ), A_new
end

function flow_μ!(α, kick, x)
    
    A = SMatrix{2,2}(1.0I)

    for _=1:kick  
        x, A = ∇H_μ!(α, x, A)
    end
    x, A
end


# dynamics returning everything we need to compute semiclassical stuff. Used in get_means_osc
function ∇H_ρ!(α, x, s, A, μ)   
    # map
    M = x[1] - α*cos(x[2])        
    Θ = mod2pi(x[2] + α*M)      
    
    # action
    s  += α*(0.5*M^2 - sin(x[2]))
    
    # monodromy
    A_loc = SMatrix{2,2}([   1.0         α*sin(x[2])      ; 
                             α      1.0 + α^2*sin(x[2])   ])
    
    A_old = copy(A)
    A     = mul!(similar(A), A_loc, A_old)
    
    # Maslov index
    i,j = 2,2
    
    μ += (A[i,j] < 0.0) & (A_old[i,j] > 0.0) ? 1 : 0
    μ += (A[i,j] > 0.0) & (A_old[i,j] < 0.0) ? 1 : 0

    SVector{2}(M, Θ), s, A, μ
end

function flow_ρ!(α, kick, x)
    s = 0.0
    A = SMatrix{2,2}(1.0I)
    μ = 0
    
    for _=1:kick  
        x, s, A, μ = ∇H_ρ!(α, x, s, A, μ)
    end
    first(x), last(x), s, A, μ
end;
