### Extracting the connected set

function findLCS(df)

    f2f = @chain df begin
        @select(:i,:j,:t)
        @orderby :i :t
        transform([:i,:j] .=> lead)
        @subset :i .== :i_lead
        @select!(:j,:j_lead)
        unique(_)
        @transform :value = 1
    end;

    f2f2 = rename(f2f,[:j_lead,:j,:value]);
    adj = unique(vcat(f2f,f2f2));
    adj = sparse(adj.j,adj.j_lead,adj.value);

    # Then create the connected set and drop firms not in the set
    conn = connected_components(SimpleGraph(adj));
    keep = conn[1];
    LCS = filter(row -> row.j in keep, df);

    return LCS
end