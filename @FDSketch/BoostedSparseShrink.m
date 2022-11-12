function Bprime = BoostedSparseShrink(obj,Aprime,l,dirac)
    [m,d] = size(Aprime);
    assert(l<=m, "l must be leq m");

%     n = 0;
%     nmax = 10;
%     while (true && n < nmax)
    while (true)
        Bprime = SparseShrink(Aprime,l);
        alpha = 6/41; % black magic
        delta = ( norm(Aprime,'Fro')^2 - norm(Bprime,'Fro')^2 ) / (alpha * l);
        C = (Aprime'*Aprime - Bprime' * Bprime) / (delta/2);
        if (VerifySpectral(obj,C))
            return
        end
        n = n + 1;
    end
end

function Bprime = SparseShrink(Aprime,l)
    [m,d] = size(Aprime);
    assert(l<=m, "l must be leq m");

    [Z,~,~] = svds(Aprime,l);
    P = Z'*Aprime;
    [~,Ltilde,V] = svds(P,l);
    if (obj.flag)
        Ltilde = diag( sqrt( diag(Ltilde).^2 - Ltilde(l,l)^2 * ones(l,1) ) );
    end
    Bprime = Ltilde * V';
end