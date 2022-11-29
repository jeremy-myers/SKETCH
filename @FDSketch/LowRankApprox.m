function [U,S,V] = LowRankApprox(obj,r)
    [U,S,V] = svd(obj.B,'econ');
    U = U(:,1:r);
    S = diag(diag(S(1:r)));
    V = V(:,1:r);
end
