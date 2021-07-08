% README
%
% Author       : Frank E. Curtis
% Description  : README file.
% Last revised : 21 June 2010

% Please cite:
%   F. E. Curtis, ``A Penalty-Interior-Point Method for Large-Scale
%   Nonlinear Optimization,'' Technical Report 10T-006, Department of
%   Industrial and Systems Engineering, Lehigh University, Bethlehem,
%   PA, USA, 2010.

% This code runs a penalty-interior-point algorithm for smooth constrained
% optimization.  The problem of interest is formulated as
%   min f(x) subject to cE(x) = 0 and cI(x) <= 0,
% where f, cE and cI are assumed to be continuously differentiable on R^n.

% The code accepts AMPL input.  A user need only provide the location of an
% AMPL .nl file in the input structure i; e.g., a typical call may be
%   >> i.nl = 'model.nl';
%   >> r = run_driver(i);
% A spamfunc.mexglx file needs to be present; please see
%   D. M. Gay, ``Hooking Your Solver to AMPL,'' Technical Report, Computing
%   Sciences Research Center, Bell Laboratories, Murray Hill, NJ, USA, 1997.

% Parameters/options can also be provided in the input structure i.  A list
% of parameters/options, along with their default values and domains, is:
%   i.algorithm    =           1; % Algorithm option (0 ~ conservative, 1 ~ aggressive)
%   i.warnings     =           0; % Warning level (0 ~ all warnings suppressed)
%   i.output_file  = 'pipal.out'; % Output file handle
%   i.verbosity    =           0; % Output verbosity level (0 ~ standard, 1 ~ detailed)
%   i.opt_err_tol  =       1e-06; % Optimality error tolerance (>0)
%   i.iter_max     =        1000; % Iteration limit (>0, integer)
%   i.grad_max     =       1e+00; % Gradient norm max, for scaling purposes (>0)
%   i.infeas_max   =       1e+02; % Infeasibility limit (>0) 
%   i.ls_factor    =       5e-01; % Line search steplength reduction factor (>0, <1)
%   i.ls_thresh    =       1e-08; % Line search sufficient decrease constant (>0, <1)
%   i.ls_frac      =       1e-02; % Line search fraction-to-boundary constant (>0, <1)
%   i.pivot_thresh =       1e-12; % Pivot threshold, for matrix factorization (>0)
%   i.slack_min    =       1e-20; % Slack variable minimum, to avoid slack=0 (>0)
%   i.shift_min    =       1e-12; % Hessian modification shift minimum value (>0)
%   i.shift_factor =       5e-01; % Hessian modification shift reduction factor (>0, <1)
%   i.shift_max    =       1e+04; % Hessian modification shift maximum value (>0)
%   i.rho_init     =       1e-01; % Penalty parameter initial value (>0)
%   i.rho_min      =       1e-12; % Penalty parameter minimum value (>=0)
%   i.rho_factor   =       5e-01; % Penalty parameter reduction factor (>0, <1)
%   i.rho_trials   =           4; % Penalty parameter number of trial values (>=1, integer)
%   i.mu_init      =       1e-01; % Interior-point parameter initial value (>0)
%   i.mu_min       =       1e-12; % Interior-point parameter minimum value (>=0)
%   i.mu_factor    =       1e-01; % Interior-point parameter reduction factor (>0, <1)
%   i.mu_trials    =          10; % Interior-point parameter number of trial values (>=1, integer)
%   i.mu_max_exp   =           0; % Interior-point parameter maximum exponent for increases (>=0, integer)
%   i.mu_max       =       1e-01; % Interior-point parameter maximum value (>0)
%   i.opt_err_mem  =           6; % Optimality error memory (>=1, integer, irrelevant if i.mu_max_exp=0)
%   i.update_con_1 =       1e-02; % Parameter updating constant 1 (>0, <1)
%   i.update_con_2 =       1e-02; % Parameter updating constant 2 (>0, <1)
%   i.update_con_3 =     1+1e-02; % Parameter updating constant 3 (>=1)

% Descriptions of the output columns are as follows:
%   Iter.        % Iteration number
%   Objective    % Objective value
%   Infeas.      % Infeasibility measure
%   Pen. Par.    % Penalty parameter value
%   I.P. Par.    % Interior-point parameter value
%   Opt. Error   % Optimality error for nonlinear optimization problem
%   Merit        % Merit function value
%   P.I.P. Error % Optimality error for penalty-interior-point subproblem 
%   Shift        % Hessian modification shift value
%   ||P.Step||   % Norm of primal search direction
%   ||D.Step||   % Norm of dual search direction
%   Lin. Red.    % Reduction in linear model of merit function
%   Quad. Red.   % Reduction in quadratic model of merit function
%   Quality      % Quality function value
%   Pri. Step.   % Primal steplength
%   Dual Step.   % Dual steplength
