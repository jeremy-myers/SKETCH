function output = RowIterator(A,model,r,k,s,l,field)
    
    [n,d] = size(A);

    %% Stream A & sketch with different models
    % Initialize models
    SketchThree = ThreeSketch(model,n,d,k,s,field);
    SketchFD_Sh = [];
    SketchFD_Ss = [];
    SketchFD_Vh = [];
    SketchFD_Vs = [];
    SketchFD_par = cell(ceil(n/l),1);
    
    % Streaming iterator
    i = 1;
    j = 0;
    while (i <= n)
        t = tic;
        if (i+l-1 <= n), idx = i:i+l-1; else, idx = i:n; end
    
        i = i+l;
        j = j+1;        
    
        % Get next rows
        Ai = A(idx,:);        
       
        %% Perform SketchThree update
        SketchThree.LinearUpdateRow(Ai,idx,1,1);
    
        %% Perform FD update
        [~, SketchFD_Sh, SketchFD_Vh] = svd( [SketchFD_Sh*SketchFD_Vh'; Ai], 'econ');
        [~, SketchFD_Ss, SketchFD_Vs] = svd( [SketchFD_Ss*SketchFD_Vs'; Ai], 'econ');
        SketchFD_Ss = diag(SketchFD_Ss);
        SketchFD_Ss = sqrt( max(SketchFD_Ss.^2 - SketchFD_Ss(l)^2,0) );
        SketchFD_Ss = diag(SketchFD_Ss);

        %% Store rows for parallel FD
        SketchFD_par{j} = Ai;

%         fprintf('Window %d of %d completed. Elapsed time: %f\n',j,ceil(n/l),toc(t));
%         semilogy(1:r,SketchFD_singVals(:,j),'b-o',1:r,SketchThree_singVals(:,j),'r-o');
%         drawnow;
%         hold on;
    end
    
    %% Merge parallel FD
    [SketchFD_Spar,SketchFD_Vpar] = merge(SketchFD_par,r,false);
    
    %% Output results
    output = struct();
    
    % Hard thresholding
    Uh = A*SketchFD_Vh;
    for i = 1:size(Uh,2), Uh(:,i) = 1/SketchFD_Sh(i,i)*Uh(:,i); end
    output.SketchFD_Uh = Uh;
    output.SketchFD_Sh = SketchFD_Sh;
    output.SketchFD_Vh = SketchFD_Vh;

    % Soft thresholding
    Us = A*SketchFD_Vs;
    for i = 1:size(Us,2)
        if (SketchFD_Ss > eps)
            Us(:,i) = 1/SketchFD_Ss(i,i)*Us(:,i);
        end        
    end
    output.SketchFD_Us = Us;
    output.SketchFD_Ss = SketchFD_Ss;
    output.SketchFD_Vs = SketchFD_Vs;

    % Parallel FD (hard thresholding)
    Upar = A*SketchFD_Vpar;
    for i = 1:size(Upar,2), Upar(:,i) = 1/SketchFD_Spar(i,i)*Upar(:,i); end
    output.SketchFD_Upar = Upar;
    output.SketchFD_Spar = SketchFD_Spar;
    output.SketchFD_Vpar = SketchFD_Vpar;

    % 3Sketch-FD hybrid
    output.SketchThree = SketchThree;  

end

function varargout = merge(C,r,trunc)
    [~,S,V] = cellfun(@(x) svd(x,'econ'), C, 'UniformOutput', false);
    B = cellfun(@(A,B) sketch(A,B,r,trunc), S, V, 'UniformOutput', false);
    [~,S,V] = svd(vertcat(B{:}), 'econ');

    switch nargout
        case 1
            varargout{1} = sketch(S,V,r,trunc);
        case 2
            [varargout{1},varargout{2}] = sketch(S,V,r,trunc);
        otherwise
            error("Too many outputs specified.");
    end
end

function varargout = sketch(S,V,r,trunc)
    switch nargout
        case 1
            if trunc
                varargout{1} = S(1:r,1:r) * V(:,1:r)';
            else
                varargout{1} = S*V';
            end
        case 2
            if trunc
                varargout{1} = S(1:r,1:r);
                varargout{2} = V(:,1:r);
            else
                varargout{1} = S;
                varargout{2} = V;
            end
        otherwise
            error("Too many outputs specified");
    end
end