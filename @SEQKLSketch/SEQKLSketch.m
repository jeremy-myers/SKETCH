classdef SEQKLSketch < matlab.mixin.SetGet

    properties
        n
        d
        k
        SEQKL_U
        SEQKL_V
        S
        ttype
        mode
        fttype
        reortho
        % ints
        verbosity 
        kmin
        kmax
        lmin
        lmax
        kstart
        fkmin
        fkmax
        extrak
        numpasses
        paramtest 
        debug 
        earlystop
        % floats
        fthresh
        % arrays
        snaps
        Utest
        Vtest
        Stest
        whch
    end
    
    %% methods
    methods
        % Constructor
        % function obj = SEQKLSketch(A, k, l, Field, Orthogonalization)
        % function obj = SEQKLSketch(A, k, l, varargin)
        % function obj = SEQKLSketch(Model, m, n, k, l, ...)
        function obj = SEQKLSketch(n,d,kmax)
            obj.n = n;
            obj.d = d;
            obj.k = 0;
            obj.kmax = kmax;

            obj.ttype = 'rel';
            obj.mode = 'restart';
            obj.fttype = [];
            obj.reortho = 'yes';
            % ints
            obj.verbosity = 0;
            obj.kmin = 1;
            obj.lmin = 1;
            obj.lmax = round(obj.kmax/sqrt(2));
            obj.kstart = obj.lmax;
            obj.fkmin = 1;
            obj.fkmax = obj.kmax;
            obj.extrak = 0;
            obj.numpasses = 1;
            obj.paramtest = 0;
            obj.debug = 0;
            obj.earlystop = obj.n;
            % floats
            obj.fthresh = [];
            % arrays
            obj.snaps = [];
            obj.Utest = [];
            obj.Vtest = [];
            obj.Stest = [];
            obj.whch = 'L';
        end
        
        % Other methods
        LinearUpdate(obj, varargin)                 % Algorithm 2 in [MR]
        % [Q, W] = SimpleLowRankApprox(obj)           % Algorithm 3
        [Q, W] = LowRankApprox(obj)                 % Algorithm 4
        % [U, S] = LowRankSymApprox(obj)              % Algorithm 5
        % [U, D] = LowRankPSDApprox(obj)              % Algorithm 6
        [Q, Sigma, V] = FixedRankApprox(obj, r)     % Algorithm 7
        % [U, D] = FixedRankSymApprox(obj, r)         % Algorithm 8
        % [U, D] = FixedRankPSDApprox(obj, r)         % Algorithm 9
        
        %% Property set methods
        % function obj = set.X(obj,value)
    %         if isequal(size(value), size(obj.X)) || isempty(obj.X)
    %             obj.X = value;
    %         else
    %             error('Size of input does not match with sketch size.')
    %         end
    %     end
        
    %     function obj = set.Y(obj,value)
    %         if isequal(size(value), size(obj.Y)) || isempty(obj.Y)
    %             obj.Y = value;
    %         else
    %             error('Size of input does not match with sketch size.')
    %         end
    %     end
        
    %     %% Methods about storage
    %     function nz = nnzTestMatrix(obj)
    %         nz = nnz(obj.Omega) + nnz(obj.Upsilon);
    %     end
        
    %     function nz = nnzSEQKLSketch(obj)
    %         nz = nnz(obj.Y) + nnz(obj.X);
    %     end
        
    %     function nz = nnz(obj)
    %         nz = nnzSEQKLSketch(obj) + nnzTestMatrix(obj);
    %     end
        
    %     function nz = numelTestMatrix(obj)
    %         nz = numel(obj.Omega) + numel(obj.Upsilon);
    %     end
        
    %     function nz = numelSEQKLSketch(obj)
    %         nz = numel(obj.Y) + numel(obj.X);
    %     end
        
    %     function nz = numel(obj)
    %         nz = numelSEQKLSketch(obj) + numelTestMatrix(obj);
    %     end
    end
end


