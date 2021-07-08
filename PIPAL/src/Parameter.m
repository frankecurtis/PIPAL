% Parameters class
classdef Parameter < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    algorithm             % Algorithm number
    opt_err_tol = 1e-06; % Optimality tolerance
    iter_max    = 1e+03; % Iteration limit
    mu_max_exp  = 0;     % Interior-point parameter maximum exponent in increases
    
  end
  
  % Class properties (constant)
  properties (Constant)
    
    grad_max     = 1e+00;   % Gradient norm limit for scaling
    infeas_max   = 1e+02;   % Infeasibility limit for penalty parameter update
    opt_err_mem  = 6;       % Optimality error history length
    ls_factor    = 5e-01;   % Line search reduction factor
    ls_thresh    = 1e-08;   % Line search threshold value
    ls_frac      = 1e-02;   % Line search fraction-to-boundary constant
    pivot_thresh = 5e-01;   % Pivot threshold for LDL factorization
    slack_min    = 1e-20;   % Slack variable bound
    shift_min    = 1e-12;   % Hessian shift (nonzero) minimum value
    shift_factor = 5e-01;   % Hessian shift update value
    shift_max    = 1e+04;   % Hessian shift maximum value
    rho_init     = 1e-01;   % Penalty parameter initial value
    rho_min      = 1e-12;   % Penalty parameter minimum value
    rho_factor   = 5e-01;   % Penalty parameter reduction factor
    rho_trials   = 4;       % Penalty parameter number of trial values per iteration
    mu_init      = 1e-01;   % Interior-point parameter initial value
    mu_min       = 1e-12;   % Interior-point parameter minimum value
    mu_factor    = 5e-01;   % Interior-point parameter reduction factor
    mu_trials    = 10;       % Interior-point parameter number of trial values per iteration
    mu_max       = 1e-01;   % Interior-point parameter maximum value
    mu_max_exp0  = 0;       % Interior-point parameter maximum exponent in increases (default)
    update_con_1 = 1e-02;   % Steering rule constant 1
    update_con_2 = 1e-02;   % Steering rule constant 2
    update_con_3 = 1e-02+1; % Adaptive interior-point rule constant
    
  end
  
  % Class methods
  methods

    % Constructor
    function p = Parameter(a)
      p.algorithm = a;
    end
    
    % Reset interior-point parameter maximum exponent in increases to default
    function resetMuMaxExp(p)
      p.mu_max_exp = p.mu_max_exp0;
    end
    
    % Set interior-point parameter maximum exponent in increases to zero
    function setMuMaxExpZero(p)
      p.mu_max_exp = 0;
    end
    
  end
  
end
