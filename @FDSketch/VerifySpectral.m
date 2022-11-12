function bool = VerifySpectral(obj,C)
    bool = true;
    obj.i = obj.i + 1;  
    dirac = obj.dirac / (2 * obj.i^2 );
    c = 2; % constant
    n = ceil(c*log10(obj.d/dirac));
    x = randn(obj.d,1);
    for i = 1:n
        x = C * x;
    end
    if norm(x) <= 1, return; else, bool = false; end
end


% ((A'*A - B'*B) / (delta/2)) * x
% = (2 / delta) * (A'*A - B'*B ) * x 
% = (2 / delta) * (A'*A*x - B'*B*x)
% = (2 / delta) * (A'*(A*x) - B'*(B*x))
% = (2 / delta) * (A'*(A*x)) - (2 / delta) * (B'*(B*x))