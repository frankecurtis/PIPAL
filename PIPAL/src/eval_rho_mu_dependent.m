function [c,z] = eval_rho_mu_dependent(c,p,q,z)

% function [c,z] = eval_rho_mu_dependent(c,p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates quantities dependent on the penalty
%                and interior-point parameters.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Evaluate slack variables
z = eval_slacks(p,q,z);

% Evaluate feasibility violations
z = eval_infeasibility(q,z);

% Evaluate merit
z = eval_merit(q,z);

% Evaluate optimality errors
z = eval_kkt_errors(q,z);

% Evaluate Hessian
[c,z] = eval_hessian(c,q,z);

% Evaluate primal-dual system matrix
[c,z] = eval_kkt_matrix(c,p,q,z);
