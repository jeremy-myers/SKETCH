function LinearUpdateRow( obj, varargin )
%LINEARUPDATE Single-View Sketch: Linear Update
%This function implements the linear updates explained in [TYUC2019]. 
%See our main reference for more details.
%
%   Require: (m x n) dimensinal update matrix (H)
%            scalars (theta) and (tau)
%   Ensure:  Modifies sketch (Y,W) to reflect linear update
%                    A = theta*A + eta*H
%            which takes the form
% 					 X = theta*X + tau*Upsilon*H
%                    Y = theta*Y + tau*H*Omega'
%                    Z = theta*Z + tau*Psi*H*phi'
%
%   S = LinearUpdate(H) updates the sketches X, Y and Z by choosing
%	(tau = 1) and (theta = 0).
%
%   S = LinearUpdate(H, tau) updates the sketches X, Y and Z by
%	choosing (theta = 1 - tau).
%
%   S = LinearUpdate(H, theta, tau) updates the sketches X, Y and Z.
%
%   S = LinearUpdate(U, V, theta, tau) where U and V are tall
%   matrices updates the sketches X, Y and Z from the factors H = U*V'
%   in an efficient way.
%
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: March 01, 2018
%Last modified: February 11, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

narginchk(2,6);

m = size(obj.Y,1); %#ok
n = size(obj.X,2);

if nargin == 3
    H   = varargin{1};
    idx = varargin{2};
    eta = 0;
    nu = 1;
    obj.Upsilon.subview = idx;
    obj.Omega.subview = idx;   
    obj.Phi.subview = idx;
    obj.Psi.subview = idx;
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y(idx,:) = eta*obj.Y(idx,:) + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
elseif nargin == 4
    H     = varargin{1};
    idx   = varargin{2};
    nu    = varargin{3};
    eta   = 1 - nu;
    obj.Upsilon.subview = idx;
    obj.Phi.subview = idx;
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y(idx,:) = eta*obj.Y(idx,:) + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
elseif nargin == 5
    H      = varargin{1};
    idx    = varargin{2};
    eta    = varargin{3};
    nu     = varargin{4};
    obj.Upsilon.subview = idx;
    obj.Phi.subview = idx;
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y(idx,:) = eta*obj.Y(idx,:) + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
else
    Hforw  = varargin{1};
    Hback  = varargin{2};
    idx    = varargin{3};
    if (~isempty(idx))
        warning("Not implemented yet.");
    end
    eta    = varargin{4};
    nu     = varargin{5};
    obj.X = eta*obj.X + nu*((obj.Upsilon*Hforw)*Hback');
    obj.Y = eta*obj.Y + nu*(Hforw*(Hback'*obj.Omega'));
    obj.Z = eta*obj.Z + nu*((obj.Phi*Hforw)*(Hback'*obj.Psi'));
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*((obj.Theta*Hforw)*Hback');
    end
end


end

