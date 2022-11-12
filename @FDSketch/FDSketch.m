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
            n = varargin{1};
            d = varargin{2};
%             k = varargin{3};
            l = varargin{3};
            if (nargin < 4)
                flag = false;
            else
                flag = varargin{4};
            end

            obj.B = zeros(l,d);
            obj.S = zeros(l,l);
            obj.V = zeros(d,l);
            obj.dirac = 0.90;
            obj.i = 0;
            obj.n = n;
            obj.d = d;
%             obj.k = k;
            obj.l = l;
            obj.flag = flag;
        end
        
        LinearUpdate(obj,varargin);
        bool = VerifySpectral(obj,varargin);
        [S,V] = LowRankApprox(obj,varargin);
    end
end