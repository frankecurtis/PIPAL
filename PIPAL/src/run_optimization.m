function r = run_optimization(i)

% function r = run_optimization(i)
%
% Author       : Frank E. Curtis
% Description  : Runs optimization algorithm; initializes
%                algorithm and runs iteration loop.
% Input        : i ~ inputs
% Output       : r ~ returns
% Last revised : 21 June 2010

% Initialize output
o = init_output(i);

% Initialize counters
c = init_counters;

% Initialize parameters
p = init_parameters(i);

% Initialize quantities
q = init_quantities(i);

% Initialize iterate
[c,z] = init_iterate(c,p,q);

% Initialize returns
r = init_returns(p,q,z);

% Print header
print_header(o,q,z);

% Print break
print_break(o,c);

% Iteration loop
while run_termination_check(r)
  
  % Print iterate information
  print_iterate(o,c,z);

  % Run search direction computation
  [c,z,d] = run_direction_computation(c,p,q,z);

  % Print search direction information
  print_direction(o,z,d);

  % Run step acceptance strategy
  [c,z,a] = run_step_acceptance(c,p,q,z,d);

  % Print step acceptance information
  print_step_acceptance(o,a);
  
  % Update iterate
  [c,z] = update_iterate(c,p,q,z,d,a,1);
  
  % Update counters
  c = update_counters(c);

  % Update returns
  r = update_returns(c,p,q,z,r);

  % Print break
  print_break(o,c);

end

% Print footer
print_footer(o,c,z,r);

% Run termination
run_termination(o);
