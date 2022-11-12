function LinearUpdate(obj, varargin)
    Aprime = varargin{1};
    l = obj.l;
    d = obj.d;
    dirac = obj.dirac;
    if (nnz(Aprime) >= (l*d))
        Bprime = BoostedSparseShrink(obj,Aprime,l,dirac);
        DenseShrink(obj,Bprime,l);
    end
end