function [U,S,V] = LowRankApprox(obj,r)
    obj.BoostedSparseShrink();
    obj.DenseShrink();
end
