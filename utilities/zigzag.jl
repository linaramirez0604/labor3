
### Estimating the worker and firm FEs using the Zig Zag method
function zigzag(df, startvals = [0,0])

    df = @chain df begin
        @select(:i,:j,:α,:ψ,:t,:lw)
        @transform :alpha_hat = startvals[1]
        @transform :psi_hat = startvals[2]
    end

    diff = Inf
    mses = fill(Inf,1,1)

    while diff > 1e-9

        # Update worker FEs
        df = @chain df begin
            @transform :lw_net = :lw - :psi_hat
            groupby(:i)
            transform(:lw_net => mean => :alpha_hat)
        end

        # Update firm FEs
        df = @chain df begin
            @transform :lw_net = :lw - :alpha_hat
            groupby(:j)
            transform(:lw_net => mean => :psi_hat)
        end

        # Calculate MSE
        mse_new = mean((df.lw .- df.alpha_hat .- df.psi_hat).^2)
        diff = mses[end] - mse_new
        mses = vcat(mses,mse_new)

    end

    return df
end