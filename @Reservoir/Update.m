function Update(obj,varargin)
    %UPDATE Summary of this function goes here
    %   Detailed explanation goes here
    A = varargin{1};
    m = size(A,1);
    if ~obj.initd
        obj.sampMat(1:obj.nsamp,:) = A(1:obj.nsamp,:);
        obj.sampids(1:obj.nsamp) = 1:obj.nsamp;
        obj.count = obj.nsamp;
        for i = obj.nsamp+1:m
            obj.count = obj.count+1;
            nu = randi(obj.count);
            if (nu <= obj.nsamp)
                obj.sampMat(nu,:) = A(i,:);
                obj.sampids(nu) = i;
            end
        end
        obj.initd = true;
    else
        for i = 1:m
            obj.count = obj.count+1;
            nu = randi(obj.count);
            if (nu <= obj.nsamp)
                obj.sampMat(nu,:) = A(i,:);
                obj.sampids(nu) = i;
            end
        end
    end
end
