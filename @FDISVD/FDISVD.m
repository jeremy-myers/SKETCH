classdef FDISVD < matlab.mixin.SetGet
    %% properties
    properties %(Access = private)
        S           % Spectral sketch
        V           % Range sketch
        m
        n
        k
        q
        svals
        sampler     % (Optional) Sampler
        error_reservoir       % (Optional) Sampler error
        fh          % (Optional) Matrix function handle
        Theta       % (Optional) Error test matrix
        W           % (Optional) Error sketch
        error_estimate % TYUC estimate
    end

    %% methods
    methods
        % Constructor
        % function obj = FDISVD(m, n, k, ...)
        function obj = FDISVD(varargin)
            m = varargin{1}; % Input size: A is (m x n)
            n = varargin{2}; % Input size: A is (m x n)
            k = varargin{3}; % range parameter
            Afun = [];
            q = 0;
            field = 'real';
            for t = 4:length(varargin)
                if isa(varargin{t},'function_handle')
                    Afun = varargin{t};
                elseif isnumeric(varargin{t})
                    q = varargin{t};
                elseif strcmp('real',varargin{t}) || strcmp('complex',varargin{t})
                    field = varargin{t};
                end
            end

            obj.S = [];
            obj.V = [];
            obj.svals = [];
            obj.error_reservoir = [];
            obj.error_estimate = [];
            obj.m = m;
            obj.n = n;
            obj.k = k;

            if ~isempty(Afun)
                obj.fh = Afun;
            end

            if q
                obj.Theta = Gauss(q,m,field);
                obj.W = zeros(q,n);
                obj.q = q;
                obj.sampler = Reservoir(q,n);
            end
        end

        % Other methods
        Update(obj, varargin);
        Stream(obj, varargin);
        [S, V] = get(obj, varargin);

        function reset(obj)
            obj.S = [];
            obj.V = [];
            obj.svals = [];
            obj.error_reservoir = [];
            obj.error_estimate = [];
            obj.Theta = Gauss(q,m,field);
            obj.W = zeros(q,n);
            obj.sampler = Reservoir(obj.q, obj.n);
        end

    end
end


