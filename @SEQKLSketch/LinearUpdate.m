function LinearUpdate(obj,varargin)

    A = varargin{1};

    if (obj.k == 0)
        obj.mode = 'expand';
    end

    if isequal(obj.mode,'expand')
        seqkl_stdpass(obj,A);
    elseif isequal(obj.mode,'restart'),
        % Perform restarted pass: A*[V ...]
        seqkl_restart(obj,A);                        ...
     elseif isequal(obj.mode,'sda'),
        % Perform impsd pass: A*[V Wgrad ...]
        seqkl_sda(obj,A);
     elseif isequal(obj.mode,'sdb'),
        % Perform impsd pass: A*[V Wgrad]
        seqkl_sdb(obj,A);
     end
end