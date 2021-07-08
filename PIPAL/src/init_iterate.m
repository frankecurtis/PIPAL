function [c,z] = init_iterate(c,p,q)

% function [c,z] = init_iterate(c,p,q)
%
% Author       : Frank E. Curtis
% Description  : Initializes iterate.  Primal iterate set by AMPL model;
%                equality constraint multipliers set to zero; inequality
%                constraint multipliers set to 1/2; evaluates initial
%                problem function values; stores initial data.
% Input        : c ~ counters
%                p ~ parameters
%                q ~ quantities
% Output       : c ~ updated counters
%                z ~ iterate
% Last revised : 21 June 2010

% Initialize penalty parameter
z.rho = p.rho_init;

% Initialize interior-point parameter
z.mu = p.mu_init;

% Initialize primal iterate
z.x = q.x_0;

% Initialize equality constraint multipliers
if q.nE > 0, z.lambdaE = zeros(q.nE,1); end;

% Initialize inequality constraint multipliers
if q.nI > 0, z.lambdaI = (1/2)*ones(q.nI,1); end;

% Evaluate first derivatives
[c,z] = eval_gradients(c,q,z,0);

% Evaluate scalings
z = eval_scalings(p,q,z);

% Scale gradients
[c,z] = eval_gradients(c,q,z,1);

% Evaluate functions
[c,z] = eval_functions(c,q,z);

% Evaluate penalty and interior-point parameter dependent quantities
[c,z] = eval_rho_mu_dependent(c,p,q,z);

% Store initial infeasibility
z.v_0 = z.v;

% Initialize last infeasibility
z.v_last = z.v;

% Initialize optimality error memory
z.kkt_last = inf*ones(p.opt_err_mem,1);

% Initialize last primal step norm
z.dx_norm_last = max([norm(z.x),norm(z.g),1]);
