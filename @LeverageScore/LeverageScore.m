classdef LeverageScore < matlab.mixin.SetGet
  %LEVERAGESCORE Summary of this class goes here
  %   Detailed explanation goes here

  properties (Access=public)
    Xi
    eta
    epsilon
    delta
  end

  properties (Access=private)
    size
  end

  methods
    function obj = LeverageScore(n, varargin)
      %RESERVOIR Construct an instance of this class
      %   Detailed explanation goes here
      obj.Xi = zeros(0,n);
      obj.eta = zeros(0,1);
      obj.epsilon = 1e-1;
      obj.delta = 1e-10;
      obj.size = 0;
    end

    Update(obj, varargin);
    [rnorms] = ResidualEstimate(obj, U, S, V, varargin);

    %% Property set methods
    function set.Xi(obj,value)
      obj.Xi = value;
    end

  end
end

