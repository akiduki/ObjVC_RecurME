function [S_dual, dt_dual, iter, Y] = twoFrm_inner_ialm(Dprime, J, lambda, tol, maxIter, grpmode)
% Inner loop of the L1 global motion estimation solver
% Code based on RASL code by Yigang Peng but solving different ADMM setup
% argmin_{S, dt} ||S||_1 s.t. J*dt + S = Dprime

% Modified with significant revision in group sparsity mode

if nargin < 3
    error('Too few arguments') ;
end

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

if nargin == 6
    % use the 6th argument to pass group sparsity parameters
    grpFlag = 1;
    % check SPAM documentation for assigning the fields in either structs
    grpGraph = grpmode.graph;
    grpParam = grpmode.param;
else
    grpFlag = 0;
end

DISPLAY_EVERY = 0;

% initialize
Y = Dprime;
% [U S V] = svd(Dprime);
% Y = U * V';
norm_two = norm(Y, 2);
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
% norm_two = norm_two / dual_norm;
% norm_inf = norm_inf / dual_norm;
dt_dual = zeros(size(J,2),1);

mu = 1/norm(Dprime) ;
rho = 1.2;

d_norm = norm(Dprime);

iter = 0;
converged = false;
while ~converged       
    iter = iter + 1;
    
    % Solve subproblem 1 w.r.t S    
    temp_T = Dprime - J*dt_dual - (1/mu)*Y;
    if grpFlag,
        % Group sparsity mode, require SPAM
        grpParam.regul = 'graph';
        grpParam.lambda = 1/mu;
        grpParam.numThreads = 4;
        S_dual = mexProximalGraph(temp_T,grpGraph,grpParam);
    else
        % Regular mode, just soft-thresholding
        S_dual = sign(temp_T) .* pos(abs(temp_T) - 1/mu);
    end
    
    % Solve subproblem 2 w.r.t dt_dual
    temp_T = Dprime - S_dual - (1/mu)*Y;
    dt_dual =  pinv(J)*temp_T;
    
    % Update Lagrangian
    Z = -Dprime + J*dt_dual + S_dual;
    Y = Y + mu*Z;
        
    mu = mu*rho;
    stoppingCriterion = norm(Z) / d_norm;
    
    if mod( iter, DISPLAY_EVERY) == 0
        disp(['#Iteration ' num2str(iter)  ...
            ' ||S||_0 ' num2str(length(find(abs(S_dual)>0)))...
            '  Stopping Criterion ' ...
            num2str(stoppingCriterion)]);
        disp(['1/mu=' num2str(1/mu)]);
    end    
    
    if stoppingCriterion <= tol
        if DISPLAY_EVERY
            disp('RASL inner loop is converged at:');
            disp(['Iteration ' num2str(iter) ...
                ' ||S||_0 ' num2str(length(find(abs(S_dual)>0)))  '  Stopping Criterion ' ...
                num2str(stoppingCriterion)]) ;
        end
        converged = true ;
    end
    
    if ~converged && iter >= maxIter
        if DISPLAY_EVERY
            disp('Maximum iterations reached') ;
        end
        converged = 1 ;       
    end
end
