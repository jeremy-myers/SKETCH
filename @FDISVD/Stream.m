function Stream(obj, varargin)
    %STREAM Summary of this function goes here
    %   Detailed explanation goes here
    A = varargin{1};
    l = varargin{2};
    m = size(A,1);
    k = obj.k;
    j = 1;   
    while (j <= m)
        if (j+l-1 <= m), idx = j:j+l-1; else, idx = j:m; end

        if isa(obj.fh,'function_handle')
            Afun = obj.fh(A(idx,:),A);
            Update(obj, Afun);
            if ~isempty(obj.sampler)
                % Reservoir sampling
                obj.sampler.Update(Afun);
                [Xi,nu1] = obj.sampler.get();
                e = zeros(k,1);
                for i = 1:k
                    e(i) = sqrt(m/obj.q)*norm(Xi*obj.V(:,i)-obj.S(i)*obj.V(nu1,i));
                end
                obj.error = [obj.error,e];

                % Random proejction sampling
                obj.W = eta*obj.W + nu1*(obj.Theta*Afun);
            end
        else
            Update(obj, A(idx,:));
            if ~isempty(obj.sampler)
                % Reservoir sampling
                obj.sampler.Update(A(idx,:));
                [Xi,nu1] = obj.sampler.get();
                e = zeros(k,1);
                for i = 1:k
                    e(i) = sqrt(m/obj.q)*norm(Xi*obj.V(:,i)-obj.S(i)*obj.V(nu1,i));
                end
                obj.error_reservoir = [obj.error_reservoir,e];

                % Random proejction sampling
                obj.W = obj.W + obj.Theta(:,idx)*A(idx,:);
                nu2 = randperm(obj.n,obj.q);
                if isreal(obj.Theta), b = 1; else, b = 2; end
                e = zeros(k,1);
                for i = 1:k
                    e(i) = (1/(b*obj.q))*norm(obj.W*obj.V(:,i)-obj.S(i)*obj.V(nu2,i));
                end
                obj.error_estimate = [obj.error_estimate,e];
            end
        end

        obj.svals = [obj.svals, obj.S];
        j = j+l;
    end
end

