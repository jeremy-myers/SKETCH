function [U,S,V] = mySketch(A,d,r)
    [m,n] = size(A);
%     Pi1 = sparse(randperm(m),1:m,sign(randn(m,1)));
%     Pi2 = sparse(randperm(m),1:m,sign(randn(m,1)));
%     Pi3 = sparse(randperm(m),1:m,sign(randn(m,1)));
%     Pi4 = sparse(randperm(m),1:m,sign(randn(m,1)));
    Pi5 = sparse(randperm(m),1:m,sign(randn(m,1)));
    Pi6 = sparse(randperm(m),1:m,sign(randn(m,1)));
    Pi7 = sparse(randperm(m),1:m,sign(randn(m,1)));
    Pi8 = sparse(randperm(m),1:m,sign(randn(m,1)));
    
    %
%     X = Pi1*A;
    X = rff(A,m,n,d)'*rff(A,m,n,d);
%     X = Pi2*X;
%     X = dct(X);
    X = rff(X,m,m,d)'*rff(X,m,m,d);
    
    %
%     Y = Pi3*A;
    Y = rff(A,m,n,d)'*rff(A,m,n,d);
%     Y = Pi4*Y;
%     Y = dct(Y);
    Y = rff(Y,m,m,d)'*rff(Y,m,m,d);
    Y = Y';
    
    %
%     Z = Pi5*A;
    Z = rff(A,m,n,d)'*rff(A,m,n,d);
%     Z = Pi6*Z;
%     Z = dct(Z);
%     Z = Pi7*Z;
%     Z = dct(Z);
%     Z = Pi8*Z;
%     Z = dct(Z);
    Z = rff(Z,m,m,d)'*rff(Z,m,m,d);
    Z = rff(Z,m,m,d)'*rff(Z,m,m,d);
    Z = rff(Z,m,m,d)'*rff(Z,m,m,d);
    Z = Z';
    
    %
    [Q,~] = qr(Y,0);
    [P,~] = qr(X',0);
    
    %
    PhiQ = Pi5*Q;
    PhiQ = dct(PhiQ);
    PhiQ = Pi6*PhiQ;
    PhiQ = dct(PhiQ);
    
    %
    PsiP = Pi7*P;
    PsiP = dct(PsiP);
    PsiP = Pi8*PsiP;
    PsiP = dct(PsiP);
    
    %
    [U1,T1] = qr(PhiQ,0);
    [U2,T2] = qr(PsiP,0);
    W = T1\(U1'*Z*U2)/T2';
%     W = T1\(U1'*Z*U1)/T1';
    
    [U,S,V] = svd(W,'econ');
    U = Q*U(:,1:r);
    S = S(1:r,1:r);
    V = P*V(:,1:r);
end

function f = rff(M,m,n,d)
      f = sqrt(2/d) * cos(randn(d,n)*M'+rand(d,1)*2*pi*ones(m,1)');
end