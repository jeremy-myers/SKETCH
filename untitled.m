clear; 

close all;
addpath(genpath('../kernels'));

%
mA = 1000; nA = 20;    % rows, cols of A
r = 5; sigmaK = 1; % rank (arbitrary), sigma of K
seedX = 0; % seed of eXperiment
kernel = GaussianKernel(sigmaK);
alpha = 1; % field is real; see TYUC19

drff_sweep = 100:100:mA;
ncol_sweep = 10:1:15;
rank_sweep = [5,10,20];
% mrkrs = {'o','x','+','^'};

for n = 1:numel(ncol_sweep)
  nA = ncol_sweep(n);

  seedX = seedX + 1;
  rng(seedX);
  A = randn(mA,nA);
  K = kernel.get(A,A);
  [V_eig,S_eig] = eig(K);
  [S_eig,I] = sort(diag(S_eig),'descend');
  S_eig = diag(S_eig);
  V_eig = V_eig(:,I);
  for k = 1:numel(rank_sweep)

    r = rank_sweep(k);
    %%
    seedX = seedX + 1;
    
    %% Truncated exact eig(K)
    S = S_eig(1:r,1:r);
    V = V_eig(:,1:r);
    
    %% TYUC19
    k3S = 4*r+alpha;   % TYUC19, Eq. 5.3
    s3S = 2*k3S+alpha;  % TYUC19, Eq. 5.3
    rng(seedX);
    Sketch_3S = ThreeSketch('SSRFT', mA, mA, k3S, s3S, 'real');
    Sketch_3S.LinearUpdate(K);
    [U_3S, S_3S, V_3S] = Sketch_3S.FixedRankApprox(r);
    
    %%
    % t = tic;
    % rng(seedX);
    % dMy = mA;
    % w = randn(dMy,nA);
    % b = rand(dMy,1)*2*pi;
    % F = sqrt(2/dMy) * cos( w*A' + b*ones(mA,1)' );
    % 
    % %
    % X = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*(F'*F);
    % X = dct(X);
    % X = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*X;
    % X = dct(X);
    % 
    % %
    % Y = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*(F'*F);
    % Y = dct(Y);
    % Y = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*Y;
    % Y = dct(Y);
    % Y = Y';
    % 
    % %
    % Z = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*(F'*F);
    % Z = dct(Z);
    % Z = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*Z;
    % Z = dct(Z);
    % Z = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*Z;
    % Z = dct(Z);
    % Z = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*Z;
    % Z = dct(Z);
    % Z = Z';
    % 
    % %
    % [Q,~] = qr(Y,0);
    % [P,~] = qr(X',0);
    % 
    % %
    % PhiQ = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*Q;
    % PhiQ = dct(PhiQ);
    % PhiQ = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*PhiQ;
    % PhiQ = dct(PhiQ);
    % PhiQ = PhiQ*sparse(randperm(mA),1:mA,sign(randn(mA,1)));
    % 
    % %
    % PsiP = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*P;
    % PsiP = dct(PsiP);
    % PsiP = sparse(randperm(mA),1:mA,sign(randn(mA,1)))*PsiP;
    % PsiP = dct(PsiP);
    % 
    % %
    % [U1,T1] = qr(PhiQ,0);
    % [U2,T2] = qr(PsiP,0);
    % W = T1\(U1'*Z*U2)/T2';
    % W = T1\(U1'*Z*U1)/T1';
    % 
    % [U_My, S_My, V_My] = svd(W,'econ');
    % U_My = Q*U_My(:,1:rK);
    % S_My = S_My(1:rK,1:rK);
    % V_My = P*V_My(:,1:rK);
    % 
    % time_My = toc(t);
    
    %%
    for d = 1:numel(drff_sweep)
      rng(seedX);
      [U_My, S_My, V_My] = mySketch(A,drff_sweep(d),r);
      norm_My(k,d)  = norm(K - U_My * (S_My*U_My')) / S_eig(1);
    end
    
    %%
%     rng(seedX);
%     w = randn(dMy,nA);
%     b = rand(dMy,1)*2*pi;
%     F = sqrt(2/dMy) * cos( w*A' + b*ones(mA,1)' );
%     
%     rng(seedX);
%     Sketch_My3S = ThreeSketch('SSRFT', mK, nK, mK, mK, 'real');
%     Sketch_My3S.LinearUpdate(F'*F);
%     [U_My3S, S_My3S, V_My3S] = Sketch_My3S.FixedRankApprox(rK);
        
    %% Compute norms       
    norm_eig(k) = norm(K- V*(S*V')) / S_eig(1);
    norm_3S(k)  = norm(K - U_3S*(S_3S*V_3S')) / S_eig(1); 
  end
  figure(nA);
  ax = gca;
  semilogy(ax, rank_sweep,norm_eig,'k-', rank_sweep, norm_3S, 'r-'); 
  hold on;
  for d = 1:numel(drff_sweep)
    semilogy(ax, rank_sweep, norm_My(:,d),'Color','b','LineStyle','-');
  end
  drawnow;
end