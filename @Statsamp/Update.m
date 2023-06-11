function Update(obj,varargin)
    %UPDATE Summary of this function goes here
    %   Detailed explanation goes here
    A = varargin{1};
    if (obj.count > 1)
        delete(obj);
    end
    insert(obj,A);
end

function insert(obj,A)
    m = size(A,1);

    Z = obj.sampMat;
    I = obj.sampids;
    s = obj.nsamp;
    x = obj.shift;
    c = obj.count;

    % Select indices and insert
    nu = randperm(m,ceil(s/c))';
    Z = [Z; A(nu,:)];
    I = [I; nu + x];
    [I,J] = sort(I,'ascend');
    Z = Z(J,:);

    % Set output
    obj.sampMat = Z;
    obj.sampids = I;
    obj.sampjds = J;

    % Increment counter & shift
    obj.count = obj.count + 1;
    obj.shift = obj.shift + m;
end

function delete(obj)
    Z = obj.sampMat;
    I = obj.sampids;
    s = obj.nsamp;
    c = obj.count;

    % Select indices and drop
    eta = randperm(s,ceil(s/c))';
    Z(eta,:) = [];
    I(eta) = [];
    [I,J] = sort(I,'ascend');
    Z = Z(J,:);

    obj.sampMat = Z;
    obj.sampids = I;
    obj.sampjds = J;
end