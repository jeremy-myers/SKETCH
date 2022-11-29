classdef FDSketch < matlab.mixin.SetGet
    properties
        B;
        S;
        V;
        dirac;
        i;
        n;
        d;
        k;
        l;
        flag;
    end
    methods
        function obj = FDSketch(varargin)
            obj.n = varargin{1};
            obj.d = varargin{2};
            obj.l = varargin{3};
            obj.flag = false;
            if (nargin == 4)
                obj.flag = varargin{4};
            end
            obj.B = zeros(obj.l,obj.d);
            obj.S = zeros(obj.l,obj.l);
            obj.V = zeros(obj.d,obj.l);
            obj.dirac = 0.90;
            obj.i = 0;
        end
        LinearUpdate(obj,varargin);
        bool = VerifySpectral(obj,varargin);
        [S,V] = LowRankApprox(obj,varargin);
    end
end