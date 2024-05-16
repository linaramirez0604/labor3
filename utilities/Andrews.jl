### Biased-corrected worker and firm variances as well as worker-firm covariance using Andrews et al. (2008)

# First let's just do the regular version before we deal with the trace thing
function andrews(df,exact=false)
    #exact: true=exact trace, false=estimated trace

    # Create the matrix form the equations
    df = @orderby(df, :i, :j, :t)
    y = df.lw
    D = indicatormat(df.i)'
    F = indicatormat(df.j)'[:,2:end] # Need to drop 1 firm dummy variable

    # Worker FEs
    M̃f = I(size(F,1)) - F*inv(F'*F)*F'
    Q̃f = M̃f
    θ̂  = inv(D' * Q̃f * D) * D' * Q̃f * y
    #df = leftjoin(df, DataFrame( i=sort(unique(df.i)), θ_hat = θ̂ ), on=:i)

    # Firm FEs
    M̃d = I(size(D,1)) - D*inv(D'*D)*D'
    Q̃d = M̃d
    ψ̂  = inv(F' * Q̃d * F) * F' * Q̃d * y
    #df = leftjoin(df, DataFrame( j=sort(unique(df.j)), ψ_hat = [0; ψ̂] ), on=:j)

    # MSE
    # mse = mean((df.lw .- df.θ_hat .- df.ψ_hat).^2)

    # Bias correction step 1: estimator for σ²ε (see FN 9 of BHLMMS)
    V = [D F]
    M = I(size(V,1)) - V * inv(V' * V) * V'
    traceM = size(df,1) - size(D,2) - size(F,2)
    σ̂²ε = y' * M * y / traceM
    
    # Infeasible bias-corrected estimators:
    Nstar = size(df,1);
    A = I(Nstar) - ones(Nstar) * inv(ones(Nstar)' * ones(Nstar)) * ones(Nstar)';
    if exact
        
        σ̂²θ = ( θ̂' * D' * A * D * θ̂ - σ̂²ε * tr( inv(D' * Q̃f * D) * D' * A * D ) ) / Nstar
        #var(df.θ_hat)

        σ̂²ψ = ( ψ̂' * F' * A * F * ψ̂ - σ̂²ε * tr( inv(F' * Q̃d * F) * F' * A * F ) ) / Nstar
        #var(df.ψ_hat)

        σ̂²θψ = ( θ̂' * D' * A * F * ψ̂ + σ̂²ε * tr( D' * A * F * inv(F' * Q̃d * F) *  F' * D * inv(D' * D) ) ) / Nstar
        #cov(df.θ_hat,df.ψ_hat)

    else
        try
            traceθ = traceApprox(D' * Q̃f * D, D' * A * D)
            σ̂²θ = ( θ̂' * D' * A * D * θ̂ - σ̂²ε * traceθ ) / Nstar
            #var(df.θ_hat)

            traceψ = traceApprox(F' * Q̃d * F, F' * A * F)
            σ̂²ψ = ( ψ̂' * F' * A * F * ψ̂ - σ̂²ε * traceψ ) / Nstar
            #var(df.ψ_hat)

            traceθψ = traceApprox2(Nstar, D'*D, F'*Q̃d*F, D'*A, F'*A, F'*A*D)    
            σ̂²θψ = ( θ̂' * D' * A * F * ψ̂ + σ̂²ε * traceθψ ) / Nstar
            #cov(df.θ_hat,df.ψ_hat) 
        catch
            (σ̂²ψ, σ̂²θ, σ̂²θψ) = (0, 0, 0)
            println("Error in calculating the Hutchingson trace estimator") 
        end
    end

    return (σ̂²ψ, σ̂²θ, σ̂²θψ)

end

function traceApprox(invmat,regmat)

    n = size(invmat,1)
    p = rand(1:n)
    
    T_p = 0
    for i = 1:p
        r_i = sign.(randn(n))
        v = invmat \ (regmat * r_i)
        T_p += r_i' * v
    end

    T_p = T_p / p

    return T_p
end


function traceApprox2(Nstar,i1,i2,r1,r2,m)
    #Adapting Gaure (2014) notation to our setting
    #Nstar is the # of worker-year observations
    #i1 is the first  inverse part
    #i2 is the second inverse part
    #r1 is the first  regular part
    #r2 is the second regular part
    #m  is the term to  compute the final estimate

    p = rand(1:Nstar)
    
    T_p = 0
    for i = 1:p
        r_i = sign.(randn(Nstar))

        v = i1 \ (r1 * r_i)
        w = i2 \ (r2 * r_i)
        
        T_p += w' * m * v
    end

    T_p = T_p / p

    return T_p
end