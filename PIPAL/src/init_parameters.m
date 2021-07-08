function p = init_parameters(i)

% function p = init_parameters(i)
%
% Author       : Frank E. Curtis
% Description  : Initializes parameters not initialized by inputs.
% Input        : i ~ inputs
% Output       : p ~ parameters
% Last revised : 21 June 2010

% Set algorithm number parameter
if isfield(i,'algorithm'),    p.algorithm    = i.algorithm;    else p.algorithm    =       1; end;

% Set termination parameters
if isfield(i,'opt_err_tol'),  p.opt_err_tol  = i.opt_err_tol;  else p.opt_err_tol  =   1e-06; end;
if isfield(i,'iter_max'),     p.iter_max     = i.iter_max;     else p.iter_max     =    1000; end;

% Set constraint/optimality tolerance parameters
if isfield(i,'grad_max'),     p.grad_max     = i.grad_max;     else p.grad_max     =   1e+00; end;
if isfield(i,'infeas_max'),   p.infeas_max   = i.infeas_max;   else p.infeas_max   =   1e+02; end;
if isfield(i,'opt_err_mem'),  p.opt_err_mem  = i.opt_err_mem;  else p.opt_err_mem  =       6; end;

% Set line search parameters
if isfield(i,'ls_factor'),    p.ls_factor    = i.ls_factor;    else p.ls_factor    =   5e-01; end;
if isfield(i,'ls_thresh'),    p.ls_thresh    = i.ls_thresh;    else p.ls_thresh    =   1e-08; end;
if isfield(i,'ls_frac'),      p.ls_frac      = i.ls_frac;      else p.ls_frac      =   1e-02; end;

% Set matrix modification parameters
if isfield(i,'pivot_thresh'), p.pivot_thresh = i.pivot_thresh; else p.pivot_thresh =   1e-12; end;
if isfield(i,'slack_min'),    p.slack_min    = i.slack_min;    else p.slack_min    =   1e-20; end;
if isfield(i,'shift_min'),    p.shift_min    = i.shift_min;    else p.shift_min    =   1e-12; end;
if isfield(i,'shift_factor'), p.shift_factor = i.shift_factor; else p.shift_factor =   5e-01; end;
if isfield(i,'shift_max'),    p.shift_max    = i.shift_max;    else p.shift_max    =   1e+04; end;

% Set penalty parameter update parameters
if isfield(i,'rho_init'),     p.rho_init     = i.rho_init;     else p.rho_init     =   1e-01; end;
if isfield(i,'rho_min'),      p.rho_min      = i.rho_min;      else p.rho_min      =   1e-12; end;
if isfield(i,'rho_factor'),   p.rho_factor   = i.rho_factor;   else p.rho_factor   =   5e-01; end;
if isfield(i,'rho_trials'),   p.rho_trials   = i.rho_trials;   else p.rho_trials   =       4; end;

% Set interior-point parameter update parameters
if isfield(i,'mu_init'),      p.mu_init      = i.mu_init;      else p.mu_init      =   1e-01; end;
if isfield(i,'mu_min'),       p.mu_min       = i.mu_min;       else p.mu_min       =   1e-12; end;
if isfield(i,'mu_factor'),    p.mu_factor    = i.mu_factor;    else p.mu_factor    =   5e-01; end;
if isfield(i,'mu_trials'),    p.mu_trials    = i.mu_trials;    else p.mu_trials    =      10; end;
if isfield(i,'mu_max_exp'),   p.mu_max_exp   = i.mu_max_exp;   else p.mu_max_exp   =       0; end;
if isfield(i,'mu_max'),       p.mu_max       = i.mu_max;       else p.mu_max       =   1e-01; end;

% Set penalty and interior-point parameter update parameters
if isfield(i,'update_con_1'), p.update_con_1 = i.update_con_1; else p.update_con_1 =   1e-02; end;
if isfield(i,'update_con_2'), p.update_con_2 = i.update_con_2; else p.update_con_2 =   1e-02; end;
if isfield(i,'update_con_3'), p.update_con_3 = i.update_con_3; else p.update_con_3 = 1+1e-02; end;
