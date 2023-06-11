classdef Reservoir < matlab.mixin.SetGet
    %RESERVOIR Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access=public)
        sampMat
        sampids
    end

    properties (Access=private)
        initd
        nsamp
        count
    end

    methods
        function obj = Reservoir(s, n, varargin)
            %RESERVOIR Construct an instance of this class
            %   Detailed explanation goes here
            obj.sampMat = [];
            obj.sampids = [];
            obj.initd = false;
            obj.nsamp = s;
            obj.count = 0;
        end

        Update(obj, varargin);
        [Xi,eta] = get(obj);
        
        function tf = isinit(obj)
            tf = obj.initd;
        end

        function sz = size(obj)
            sz = obj.nsamp;
        end

        function c = curr_row_idx(obj)
            c = obj.count;
        end

        function reset(obj)
            obj.sampMat = [];
            obj.sampids = [];
            obj.initd = false;
            obj.count = 0;
        end

    end
end
