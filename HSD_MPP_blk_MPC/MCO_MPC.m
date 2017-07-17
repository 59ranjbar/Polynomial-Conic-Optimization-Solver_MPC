function [x_opt,y_opt,s_opt] = MCO_MPC(A0,b0,c0,blk,scale_data,file_name)
% This script is an implementation of Non-Symmeteric 
% Homogeneous Self-Dual Mehrotra Predictor-Corrector
% Interior Point Method  
% for Multi-Block Moment Cone and its dual problem, this is,
% Univariate Non-Negative Polynomial 
% in Chebyshev Basis on interval  [-1 1].
% This is, min      c_1'*x_1+...+c_k'*x_k
%             s.t.      A_1*x_1+...+A_k*x_k = b,
%                               x_i  in  M^{n_i}_[-1 1]_Ch    i=1,...,k.
% where dual problem is 
%
%         max    b'*y
%          s.t.    A'*y + s_i = c_i,                       i=1,...,k
%                             s_i  in  P^{n_i^+}_[-1 1]_Ch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Details
%   A = [A_1,...,A_k]
%   c = [c_1,...,c_k]
%   b the right hand side vector
%   scale_data = 0 if scaling is not needed and 1 otherwise
%   file_name = file's name to write the iterations information
%   blk = a vector of blocks' dimensions
%%%%%%%%%%%%%%%%%%%%%%%%% 
%Initialize Pars
%%%%%%%%%%%%%%%%%%%%%%%%%
expon = 3;
gap_tol = 1e-6;              % duality gap tolerance
inf_tol = 1e-6;                % infeasibility tolerance
maxIter = 100;              % maximum number of iterations
iter = 0;                        % iteration counter for the main while loop
Stop = 0;                     % Stopping flag for main loop

if nargin < 5
    scale_data = 1 ;
end

if nargin < 6
    file_name = 'iter_MCO_MPC';
end
%%%%%%%%%%%%%%%%%%%%%%
% Initial strictly feasible points for cones, nu, mu
%%%%%%%%%%%%%%%%%%%%%%
[ x0, y0, s0, tau0, kappa0, mu0, nu ] = initial( blk, b0 );

%%%%%%%%%%%%%%%%%%%%
% Scaling A, b, c
%%%%%%%%%%%%%%%%%%%%%
if scale_data == 1
    [A,b,c,normA,normc,normb] = scaling(blk,A0,b0,c0);  
else
     A = A0; b = b0; c = c0; normb = 1; normc = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Initial Residuals
%%%%%%%%%%%%%%%%%%%%%%%
x = x0; y = y0; s = s0; kappa = kappa0; tau = tau0; mu = mu0;
rp = tau*b - A*x;
rd = tau*c - A'*y - s;  
rg = kappa + c'*x - b'*y;

%%%%%%%%%%%%%%%%%%%%%%%
%Print date and the Initial information 
%%%%%%%%%%%%%%%%%%%%%%% 
file_id = fopen(file_name , 'a');                    % where to write the iterations
fprintf(file_id, '|Ax-b|= %.2e, |A''y+s-c|= %.2e, |b''y-c''x|= %.2e, x''s= %.2e \n',...
    norm(A0*x0-b0), norm(A0'*y0+s0-c0), abs(b0'*y0-c0'*x0), x0'*s0);
fprintf(file_id, ' k:   rp            rd           |cx-by|   |x''*s|      mu            alpha_a     alpha_c      tau        kappa \n');

%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%
while(~Stop && iter < maxIter) 
	iter = iter +1;
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Affine Phase 
    %%%%%%%%%%%%%%%%%%%%%%%%
    [ dx_a, dy_a, ds_a, dtau_a, dkappa_a, schur, Hx_inv, gx, Hx_inv_At ] = ...
       affinedir( blk, x, s, tau, kappa, A, b, c, rp, rd, rg, mu);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Step-Length in Affine Phase 
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [alpha_a_pos(iter), alphaM_a, alphaP_a] = step_length( blk, ...
         x, s, tau, kappa, dx_a, ds_a, dtau_a, dkappa_a );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find step-length that z_alpha is in Neib
    %%%%%%%%%%%%%%%%%%%%%%%%%%     
    j=0; al=[]; inf_norm_al_a=[]; hess_norm_al_a=[]; mu_al_a=[];
	for al = alpha_a_pos(iter):-.1:0
        j = j+1;
%         [inf_norm_al_a(j), hess_norm_al_a(j), mu_al_a(j)] = centrality(blk,...
%              x+al*dx_a, s+al*ds_a, tau+al*dtau_a, kappa+al*dkappa_a, nu); 
        [inf_norm_al_a(j), mu_al_a(j)] = centrality2(blk,...
             x+al*dx_a, s+al*ds_a, tau+al*dtau_a, kappa+al*dkappa_a, nu); 
	end

    if isempty(find(inf_norm_al_a <= .95*mu_al_a,1))
        alpha_a(iter) = 0;
    else
        alpha_a(iter) = alpha_a_pos(iter) - .1*(find(inf_norm_al_a <= .95*mu_al_a,1)-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    % centering parameters
    %%%%%%%%%%%%%%%%%%%%%%%
    sigma(iter) = (1 - alpha_a(iter))^expon;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combined Phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [dx_c, dy_c, ds_c, dtau_c, dkappa_c] = combinedir(blk, x, s, tau, kappa, ... 
                  A, b, c, rp, rd, rg, mu, sigma(iter), schur, Hx_inv, gx, Hx_inv_At);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Step Length in Combine Phase 
    %%%%%%%%%%%%%%%%%%%%%%%%% 
    [alpha_c_pos(iter), alphaM_c, alphaP_c] = step_length( blk, ...
         x, s, tau, kappa, dx_c, ds_c, dtau_c, dkappa_c );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find step-length that z_alpha is in Neib
    %%%%%%%%%%%%%%%%%%%%%%%%%%
	al=[]; j=0; inf_norm_al_c=[]; hess_norm_al_c=[]; mu_al_c=[];
    for al = alpha_c_pos(iter):-.1:0
        j = j+1;
%         [inf_norm_al_c(j), hess_norm_al_c(j), mu_al_c(j)] = centrality(blk,...
%                  x+al*dx_c, s+al*ds_c, tau+al*dtau_c, kappa+al*dkappa_c, nu); 
        [inf_norm_al_c(j), mu_al_c(j)] = centrality2(blk,...
                 x+al*dx_c, s+al*ds_c, tau+al*dtau_c, kappa+al*dkappa_c, nu);  
    end

    alpha_c(iter) = alpha_c_pos(iter) - .1*(find(inf_norm_al_c <= .95*mu_al_c,1)-1);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Prepare points for next iteration
    %%%%%%%%%%%%%%%%%%%%%%
    x = x + alpha_c(iter)*dx_c;
    y = y + alpha_c(iter)*dy_c;
    s = s + alpha_c(iter)*ds_c;
    tau = tau + alpha_c(iter)*dtau_c;
    kappa = kappa + alpha_c(iter)*dkappa_c;
    
    rp = tau*b - A*x;
    rd =  tau*c - A'*y - s;
    rg = kappa + c'*x - b'*y;

    mu = (x'*s + tau*kappa) / nu;
        
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Chech stopping criteria 
    %%%%%%%%%%%%%%%%%%%%%%%%
    gap(iter) = (normb*normc/tau^2)*(x'*s); 
    obj(:,iter) = (normb*normc/tau)*[c'*x ;  b'*y]; 
    obj_gap = obj(1,iter) - obj(2,iter);
    
    rel_gap = gap(iter)/max(1,mean(abs(obj(:,iter))));

    pri_inf  = norm(rp)/tau; 
    dual_inf = norm(rd)/tau; 
    inf_meas(iter) = max(pri_inf, dual_inf); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Print Iteration's Information
    %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(file_id, '%2d:  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  \n',...
                        iter,  pri_inf, dual_inf, rel_gap, abs(gap(iter)),...
                        mu/tau^2, alpha_a(iter), alpha_c(iter), tau, kappa ); 

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Check the Stopping Criteria
    %%%%%%%%%%%%%%%%%%%%%%%%      
    if  max(rel_gap, inf_meas(iter)) <= 1e2*gap_tol
        Stop = 1;
         if scale_data ==1 
             for i = 1:length(blk)
                 ind_i = sum(blk(1:i-1))+1:sum(blk(1:i));

                 x_opt(ind_i,1) = x(ind_i)*(normb/(normA(i)*tau));
                 s_opt(ind_i,1) = s(ind_i)*(normc*normA(i)/tau);  
             end
             y_opt = y*(normc/(tau)); 
         else 
             x_opt = x/tau;
             y_opt = y/tau; 
             s_opt = s/tau;
         end
         fprintf(file_id, 'optimal points has been reached \n');
         fprintf(file_id, '|Ax-b|= %.2e, |A''y+s -c|= %.2e, |b''y-c''x|= %.2e, x''s= %.2e, mu= %.2e \n',...
            norm(A0*x_opt-b0), norm(A0'*y_opt+s_opt-c0), abs(b0'*y_opt-c0'*x_opt), x_opt'*s_opt, mu/tau^2);
    elseif max(mu / mu0 , (tau / kappa) / (tau0 /kappa0 )) <1e3*eps
        Stop = -1 ;
        fprintf(file_id, 'Infeasibility has been detected \n')
    end    

    
end

fprintf(file_id, '------------------------------------------------- \n')
fclose(file_id);


end









