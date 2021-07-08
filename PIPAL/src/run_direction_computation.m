function [c,z,d] = run_direction_computation(c,p,q,z)

% function [c,z,d] = run_direction_computation(c,p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Runs search direction computation; updates
%                penalty and interior-point parameters.
% Input        : c ~ counters
%                p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : c ~ updated counters
%                z ~ updated iterate
%                d ~ direction
% Last revised : 21 June 2010

% Check for interior-point parameter update based on optimality error
while z.mu > p.mu_min & z.kkt(3) <= max([z.mu;p.opt_err_tol-z.mu])
  
  % Restrict interior-point parameter increase
  p.mu_max_exp = 0;
  
  % Update interior-point parameter
  if z.mu > p.mu_min

    % Decrease interior-point
    z.mu = max(p.mu_min,p.mu_factor*z.mu);

    % Evaluate penalty and interior-point parameter dependent quantities
    [c,z] = eval_rho_mu_dependent(c,p,q,z);

  end
  
end
  
% Check for penalty parameter update based on optimality error
if (z.kkt(2) <= p.opt_err_tol & z.v > p.opt_err_tol*max(1,z.v_0)) | z.v > max([z.v_0;z.v_last;p.infeas_max])
 
  % Update penalty parameter
  if z.rho > p.rho_min

    % Decrease penalty parameter
    z.rho = max(p.rho_min,p.rho_factor*z.rho);

    % Evaluate penalty and interior-point parameter dependent quantities
    [c,z] = eval_rho_mu_dependent(c,p,q,z);

  end

end

% Set last penalty parameter
z.rho_last = z.rho;

% Check for conservative algorithm
if p.algorithm == 0

  % Evaluate primal-dual right-hand-side
  z = eval_kkt_rhs(q,z);
  
  % Evaluate direction
  d = run_kkt_solver(p,q,z);
  
  % Evaluate models
  d = eval_models(q,z,d);
  
  % Return
  return;

end

% Check KKT memory for potential mu increase limit
if z.kkt(2) > max(z.kkt_last)

  % Restrict mu increase
  p.mu_max_exp = 0;

end

% Store current penalty and interior-point parameters
rho_curr = z.rho; mu_curr = z.mu;
  
% Evaluate direction for current penalty and interior-point parameters
z.rho = rho_curr;
z.mu  = mu_curr;
z     = eval_kkt_rhs(q,z);
d11   = run_kkt_solver(p,q,z);

% Evaluate direction for zero interior-point parameter
z.rho = rho_curr;
z.mu  = 0;
z     = eval_kkt_rhs(q,z);
d10   = run_kkt_solver(p,q,z);

% Evaluate direction for zero penalty parameter
z.rho = 0;
z.mu  = mu_curr;
z     = eval_kkt_rhs(q,z);
d01   = run_kkt_solver(p,q,z);

% Set trial interior-point parameter values
Mu = max(p.mu_min,min(p.mu_factor.^([p.mu_trials-1:-1:0]-p.mu_max_exp)*mu_curr,p.mu_max));

% Initialize feasibility direction data
lred0_0_mu = zeros(1,p.mu_trials);

% Loop through interior-point parameter values
for j = 1:p.mu_trials

  % Set penalty and interior-point parameters
  z.rho = 0; z.mu = Mu(j);
    
  % Evaluate direction
  d = eval_direction(q,d11,d10,d01,(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr),z);
  
  % Cut length
  d.x = min(z.dx_norm_last/max(d.x_norm,1),1)*d.x;
    
  % Run fraction-to-boundary
  a = run_fraction_to_boundary(p,q,z,d);
    
  % Cut length
  d = run_cut_step(q,d,a);

  % Evaluate models
  d = eval_models(q,z,d);
    
  % Set feasibility direction data
  lred0_0_mu(j) = d.lred0;
      
end

% Initialize updating data
ltred0_rho_mu = zeros(p.mu_trials);
qtred_rho_mu  = zeros(p.mu_trials);
m_rho_mu      = zeros(p.mu_trials);

% Initialize check
check = 0;

% Loop through penalty parameter values
for i = 1:p.rho_trials
  
  % Set penalty parameter
  z.rho = max(p.rho_min,(p.rho_factor^(i-1))*rho_curr);

  % Set last penalty parameter
  if rho_curr > z.kkt(1)^2, z.rho_last = z.rho; end;
  
  % Loop through interior-point parameter values
  for j = 1:p.mu_trials
    
    % Set interior-point parameter
    z.mu = Mu(j);

    % Evaluate direction
    d = eval_direction(q,d11,d10,d01,(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr),z);

    % Run fraction-to-boundary
    a = run_fraction_to_boundary(p,q,z,d);
      
    % Cut steps
    d = run_cut_step(q,d,a);
      
    % Evaluate models
    d = eval_models(q,z,d);
    
    % Set updating data
    ltred0_rho_mu(j) = d.ltred0;
    qtred_rho_mu(j)  = d.qtred;
    m_rho_mu(j)      = d.m;
      
    % Check updating conditions for infeasible points
    if z.v > p.opt_err_tol*max(1,z.v_0) & (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) | qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) | z.rho > z.kkt(1)^2), m_rho_mu(j) = inf; end;

    % Check updating conditions for feasible points
    if z.v <= p.opt_err_tol*max(1,z.v_0) & qtred_rho_mu(j) < 0, m_rho_mu(j) = inf; end;

  end

  % Find minimum m for current rho
  z.m_min = min(m_rho_mu);
    
  % Check for finite minimum
  if z.m_min < inf

    % Loop through mu values
    for j = 1:p.mu_trials
    
      % Set mu
      mu = Mu(j);
      
      % Check condition
      if m_rho_mu(j) <= p.update_con_3*z.m_min, z.mu = mu; end;
        
    end
    
    % Set condition check
    check = 1;
      
    % Break loop
    break;
      
  end

end

% Check conditions
if check == 0, z.rho = rho_curr; z.mu = mu_curr; end;

% Evaluate merit
z = eval_merit(q,z);

% Evaluate primal-dual right-hand-side
z = eval_kkt_rhs(q,z);
  
% Evaluate direction
d = run_kkt_solver(p,q,z);
  
% Evaluate models
d = eval_models(q,z,d);
