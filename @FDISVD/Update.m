function Update( obj, varargin )

    A = varargin{1};
    m = size(A,1);
    k = obj.k;
    if issparse(A)
        [~,S,V] = svds([diag(obj.S) * obj.V'; A], k, 'largest', 'Tolerance', 1);
    else
        [~,S,V] = svd([diag(obj.S) * obj.V'; A]);
    end

    % Set output
    obj.S = diag(S(1:k,1:k));
    obj.V = V(:,1:k);

end
