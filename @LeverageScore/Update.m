function Update(obj, varargin)
%RIDGELEVERAGESCORE Summary of this function goes here
%   Detailed explanation goes here
  A = varargin{1};
  [n,d] = size(A);

  epsilon = obj.epsilon;
  delta = obj.delta;
  if (nargin == 4)
    epsilon = varargin{2};
    delta = varargin{3};
  end

  lambda = delta/epsilon;
  c = 8 * log(d/epsilon^2);

  I = eye(d);
  for i = 1:n
    ai = A(i,:);
    li = min((1+epsilon)*ai*((obj.Xi'*obj.Xi + lambda*I)\ai'), 1);
    pi = min(c*li, 1);
    eta = rand;
    if (pi > eta)
      obj.Xi = [obj.Xi; 1/sqrt(pi)*ai];
      obj.size = obj.size + 1;
    end
  end
end

