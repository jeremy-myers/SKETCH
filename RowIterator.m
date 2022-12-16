function output = RowIterator(A,model,r,k,s,l,field)
    
    [n,d] = size(A);

    %% Stream A & sketch with different models
    % Initialize models
    SketchThree = ThreeSketch(model,n,d,k,s,field);
    SketchFD_Sh = [];
    SketchFD_Ss = [];
    SketchFD_Vh = [];
    SketchFD_Vs = [];
    
    % Model errors
%     SketchThree_err_lra = zeros(ceil(n/l),1);
%     SketchFD_err_lra = zeros(ceil(n/l),1);

    % Singular values
%     SketchThree_singVal = zeros(r,ceil(n/l));
%     SketchFD_singVal = zeros(r,ceil(n/l));

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
%         [SketchThree_U, SketchThree_S, SketchThree_V] = SketchThree.FixedRankApprox(r);
    
        %% Perform FD update
        [~, SketchFD_Sh, SketchFD_Vh] = svd( [SketchFD_Sh*SketchFD_Vh'; Ai], 'econ');
        [~, SketchFD_Ss, SketchFD_Vs] = svd( [SketchFD_Ss*SketchFD_Vs'; Ai], 'econ');
        SketchFD_Ss = diag(SketchFD_Ss);
        SketchFD_Ss = sqrt( max(SketchFD_Ss.^2 - SketchFD_Ss(l)^2,0) );
        SketchFD_Ss = diag(SketchFD_Ss);
%         SketchFD_S = SketchFD_S(1:r,1:r);
%         SketchFD_V = SketchFD_V(:,1:r);
               
%         %% Compute error
%         SketchThree_err_lra(j) = norm(A-A*(SketchThree_V*SketchThree_V'));
%         SketchFD_err_lra(j) = norm(A-A*(SketchFD_V*SketchFD_V'));
% 
%         %% Store singular values
%         SketchThree_singVal(:,j) = diag(SketchThree_S);
%         SketchFD_singVal(:,j) = diag(SketchFD_S);

%         fprintf('Window %d of %d completed. Elapsed time: %f\n',j,ceil(n/l),toc(t));
%         semilogy(1:r,SketchFD_singVals(:,j),'b-o',1:r,SketchThree_singVals(:,j),'r-o');
%         drawnow;
%         hold on;
    end
    
    %% Output results
    output = struct();
    
    Uh = A*SketchFD_Vh;
    for i = 1:size(Uh,2), Uh(:,i) = 1/SketchFD_Sh(i,i)*Uh(:,i); end
    output.SketchFD_Uh = Uh;
    output.SketchFD_Sh = SketchFD_Sh;
    output.SketchFD_Vh = SketchFD_Vh;

    Us = A*SketchFD_Vs;
    for i = 1:size(Us,2)
        if (SketchFD_Ss > eps)
            Us(:,i) = 1/SketchFD_Ss(i,i)*Us(:,i);
        end        
    end
    output.SketchFD_Us = Us;
    output.SketchFD_Ss = SketchFD_Ss;
    output.SketchFD_Vs = SketchFD_Vs;
%     output.SketchFD_err_lra = SketchFD_err_lra(1:j);
%     output.SketchFD_singVal = SketchFD_singVal(:,1:j);  

    output.SketchThree = SketchThree;
%     output.SketchThree_U = SketchThree_U;
%     output.SketchThree_S = SketchThree_S;
%     output.SketchThree_V = SketchThree_V;
%     output.SketchThree_err_lra = SketchThree_err_lra(1:j);    
%     output.SketchThree_singVal = SketchThree_singVal(:,1:j);
  

end