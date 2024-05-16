
##### Constructing Employer-Employee matched data
###  Create a mobility matrix
function mobilityMatrix(nk, nl, csort=0.5, cnetw=0.2, csig=0.5, α_sd=1, ψ_sd=1)
    # nk    number of firm types
    # nl    number of worker types
    # csort sorting effect
    # cnetw network effect
    # csig  cross-sectional standard deviation

    # approximate each distribution with some points of support
    ψ = quantile.(Normal(), (1:nk) / (nk + 1)) * ψ_sd
    α = quantile.(Normal(), (1:nl) / (nl + 1)) * α_sd

    # Let's create type-specific transition matrices
    # We are going to use joint normals centered on different values
    G = zeros(nl, nk, nk)
    for l in 1:nl, k in 1:nk
        G[l, k, :] = pdf.( Normal(0, csig), ψ .- cnetw * ψ[k] .- csort * α[l])
        G[l, k, :] = G[l, k, :] ./ sum(G[l, k, :])
    end

    # We then solve for the stationary distribution over psis for each alpha value
    # We apply a crude fixed point approach
    H = ones(nl, nk) ./ nk
    for l in 1:nl
        M = transpose(G[l, :, :])
        for i in 1:100
            H[l, :] = M * H[l, :]
        end
    end

    return α, ψ, G, H
end

### Simulate a panel
# The next step is to simulate our network given our transition rules
function balancedPanel(λ,nt,ni,G,H)
    # λ  assuming moving probability is fixed
    # nt number of time periods
    # ni number of workers

    # We simulate a balanced panel
    ll = zeros(Int64, ni, nt) # Worker type
    kk = zeros(Int64, ni, nt) # Firm type
    spellcount = zeros(Int64, ni, nt) # Employment spell

    for i in 1:ni
        
        # We draw the worker type
        l = rand(1:nl)
        ll[i,:] .= l
        
        # At time 1, we draw from H
        kk[i,1] = sample(1:nk, Weights(H[l, :]))
        
        for t in 2:nt
            if rand() < λ
                kk[i,t] = sample(1:nk, Weights(G[l, kk[i,t-1], :]))
                spellcount[i,t] = spellcount[i,t-1] + 1
            else
                kk[i,t] = kk[i,t-1]
                spellcount[i,t] = spellcount[i,t-1]
            end
        end
    end

    return ll, kk, spellcount
end

### Attach firm ids to types
# The final step is to assign identities to the firms. We are going to do this is a relatively simple way, by simply randomly assigning firm ids to spells.
function firmIDToType(ni, nt, kk, spellcount, firms_per_type)
    jj = zeros(Int64, ni, nt) # Firm identifiers

    draw_firm_from_type(k) = sample(1:firms_per_type) + (k - 1) * firms_per_type

    for i in 1:ni
        
        # extract firm type
        k = kk[i,1]
        
        # We draw the firm (one of firms_per_type in given group)
        jj[i,1] = draw_firm_from_type(k)
        
        for t in 2:nt
            if spellcount[i,t] == spellcount[i,t-1]
                # We keep the firm the same
                jj[i,t] = jj[i,t-1]
            else
                # We draw a new firm
                k = kk[i,t]
                
                new_j = draw_firm_from_type(k)            
                # Make sure the new firm is actually new
                while new_j == jj[i,t-1]
                    new_j = draw_firm_from_type(k)
                end
                
                jj[i,t] = new_j
            end
        end
    end

    # Make sure firm ids are contiguous
    contiguous_ids = Dict( unique(jj) .=> 1:length(unique(jj))  )
    jj .= getindex.(Ref(contiguous_ids),jj);

    return jj
end


### Combine all three functions into 1
function createData(nt,              # number of time periods
                    ni,              # number of individuals
                    firms_per_type,  # number of firms per type
                    nl,              # number of worker types
                    nk,              # number of firm types
                    λ,               # probability of moving
                    csort=0.5,       # sorting effect
                    cnetw=0.2,       # network effect
                    csig=0.5,        # cross-sectional standard deviation
                    α_sd=1,          # worker effect sd
                    ψ_sd=1,          # firm effect sd
                    wages = false,   # boolean to add wages or not
                    w_sigma = 0.2)   # sd of wages

    (α, ψ, G, H) = mobilityMatrix(nk,nl,csort,cnetw,csig,α_sd,ψ_sd);
    (ll, kk, spellcount) = balancedPanel(λ,nt,ni,G,H);
    jj = firmIDToType(ni, nt, kk, spellcount, firms_per_type);

    ii = repeat(1:ni,1,nt);
    tt = repeat((1:nt)',ni,1);
    df = DataFrame(i=ii[:], j=jj[:], l=ll[:], k=kk[:], α=α[ll[:]], ψ=ψ[kk[:]], t=tt[:], spell=spellcount[:]);

    if wages
        df[!, :lw] = df.α + df.ψ + w_sigma * rand(Normal(), size(df)[1]);
    end

    return df
end