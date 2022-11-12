function DenseShrink(obj,A,l)
%     l = obj.l;
%     if (obj.l > obj.k)
%         l = obj.k;
%     end
    [m,d] = size(A);
    assert(l<=m, "l must be leq m");    
    [~,Ltilde,V] = svds([obj.B; A],l);
    if (obj.flag)
        Ltilde = diag( sqrt( diag(Ltilde).^2 - Ltilde(l,l)^2 * ones(l,1) ) );
    end
    obj.B = Ltilde * V';
end