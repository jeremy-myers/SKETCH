classdef KernelSSRFT < DimRedux
    %SSRFT implements a subclass of DimRedux for the sketching method
    %described [TYUC2019]. See our paper and its supplementary materials 
    %for more details.
    %
    %[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
    %Streaming Low-Rank Matrix Approximation with an Application to
    %Scientific Simulation. 
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Last modified: July 22, 2018
    %
    %SKETCHv1.1
	%Copyright (C) 2018 Laboratory for Information and Inference Systems
	%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	%This code is a part of SKETCH toolbox. 
	%Please read COPYRIGHT before using this file.    
    
    %% properties
    properties (Access = private)
        field
        Pi1
        Pi2
        coords
        sigma
    end
    
    %% methods
    methods
        % Constructor
        function obj = KernelSSRFT(k, n, varargin)
            obj@DimRedux(k,n);
%             obj.coords = randperm(n,k);
            obj.sigma = 1;
            if length(varargin) >= 1
                obj.field = lower(varargin{1});
            else
                obj.field = 'real';
            end
            if strcmp(obj.field,'real')
                obj.Pi1 = sparse(randperm(n),1:n,sign(randn(n,1)));
                obj.Pi2 = sparse(randperm(n),1:n,sign(randn(n,1)));
            elseif strcmp(obj.field,'complex')
                obj.Pi1 = sparse(randperm(n),1:n,...
                    sign(randn(n,1) + 1i*randn(n,1)));
                obj.Pi2 = sparse(randperm(n),1:n,...
                    sign(randn(n,1) + 1i*randn(n,1)));
            else
                err('Input ''field'' should be ''real'' or ''complex''.')
            end
        end
        
        %% Other methods
        function B = LeftApply(obj,M)
            [m,n] = size(M);
            k = size(obj,1);
            w = randn(m,k);
            b = rand(k,n)*2*pi;
            B = obj.Pi1*M;
            B = sqrt(2/k) * cos( w'*B + b );          
        end
        
        function B = RightApply(obj,M)
            [m,n] = size(X);
            k = size(obj,1);
            w = randn(k,n);
            b = rand(k,m)*2*pi;
            B = M*obj.Pi1';
            B = sqrt(2/k) * cos( B*w' + b' );

%             B = zeros(k,d);
%             B(:,obj.coords) = M;
%             if isreal(obj)
%                 B = idct(B')';
%             else
%                 B = d*ifft(B')';
%             end
%             B = B*obj.Pi2;
%             if isreal(obj)
%                 B = idct(B')';
%             else
%                 B = d*ifft(B')';
%             end
%             B = B*obj.Pi1;
        end
        
        %% Overloaded methods
        
        function flag = isreal(obj)
            flag = strcmp(obj.field,'real');
        end
        
        function flag = issparse(obj) %#ok
            flag = false;
        end
        
        function nz = nnz(obj)
            nz = numel(obj);
        end
    end
    
end

