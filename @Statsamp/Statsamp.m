classdef Statsamp < matlab.mixin.SetGet
    %RESERVOIR Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access=public)
        sampMat
        sampids
        sampjds
    end

    properties (Access=private)
        nsamp
        count
        shift
    end

    methods
        function obj = Statsamp(s, n, varargin)
            %RESERVOIR Construct an instance of this class
            %   Detailed explanation goes here
            obj.sampMat = zeros(0,n);
            obj.sampids = zeros(0,1);
            obj.sampjds = zeros(0,1);
            obj.nsamp = s;
            obj.count = 1;
            obj.shift = 0;
        end

        Update(obj, varargin);
        [Xi,nu] = get(obj);

        %% Property set methods
        function set.sampMat(obj,value)
            obj.sampMat = value;
        end

        function tf = isinit(obj)
            tf = obj.initd;
        end

        function sz = size(obj)
            sz = obj.nsamp;
        end

        function c = curr_row_idx(obj)
            c = obj.count;
        end

    end
end
