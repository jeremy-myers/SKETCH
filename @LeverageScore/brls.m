function [Ahat] = RidgeLeverageScore(obj, varargin)
%RIDGELEVERAGESCORE Summary of this function goes here
%   Detailed explanation goes here
  A = varargin{1};
  [n,d] = size(A);

  epsilon = 1e-1;
  delta = 1e-10;
  if (nargin == 4)
    epsilon = varargin{2};
    delta = varargin{3};
  end

  lambda = delta/epsilon;
  c = 8 * log(d/epsilon^2);

  if ~obj.isInitialized
    obj.S = zeros(0,d);
    obj.isInitialized = true;
  end

  I = eye(d);
%   for i = 1:n
%     ai = A(i,:);
    Shat = inv(obj.S'*obj.S + lambda*I);
    li = min((1+epsilon)*A*(Shat*A'), 1);
%     pi = min(c*li, 1);
%     eta = rand;
%     if (pi > eta)
%       obj.S = [obj.S; 1/sqrt(pi)*ai];
%     end
%   end
end

