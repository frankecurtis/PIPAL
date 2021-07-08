function [c,z,d] = update_iterate(c,p,q,z,d,a,opt)

% function [c,z,d] = update_iterate(c,p,q,z,d,a,opt)
%
% Author       : Frank E. Curtis
% Description  : Updates iterate; recomputes problem functions.
% Input        : c   ~ counters
%                p   ~ parameters
%                q   ~ quantities
%                z   ~ iterate
%                d   ~ direction
%                a   ~ step acceptance values
%                opt ~ update option
%                        0 ~ update line search quantities
%                        1 ~ update remaining quantities
% Output       : c   ~ updated counters
%                z   ~ updated iterate
%                d   ~ updated direction
% Last revised : 21 June 2010

% Update point
z = update_point(q,z,d,a);

% Evaluate line search quantities
if opt == 0
  
  % Evaluate functions
  [c,z] = eval_functions(c,q,z);

  % Evaluate slack variables
  z = eval_slacks(p,q,z);

  % Evaluate merit
  z = eval_merit(q,z);
  
% Evaluate remaining quantities
else
  
  % Update last infeasibility
  z.v_last = z.v;

  % Update last primal direction norm
  z.dx_norm_last = a.p*d.x_norm;

  % Update last primal steplength
  z.cut_last = (a.p < a.p0);

  % Evaluate first derivatives
  [c,z] = eval_gradients(c,q,z,2);

  % Evaluate penalty and interior-point parameter dependent quantities
  [c,z] = eval_rho_mu_dependent(c,p,q,z);

  % Update optimality error memory
  z.kkt_last = [z.kkt(2); z.kkt_last(1:p.opt_err_mem-1)];

end
