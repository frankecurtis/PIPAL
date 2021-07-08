function [c,z,a] = run_step_acceptance(c,p,q,z,d)

% function [c,z,a] = run_step_acceptance(c,p,q,z,d)
%
% Author       : Frank E. Curtis
% Description  : Runs step acceptance strategy.
% Input        : c ~ counters
%                p ~ parameters
%                q ~ quantities
%                z ~ iterate
%                d ~ direction
% Output       : c ~ updated counters
%                a ~ step acceptance values
% Last revised : 21 June 2010

% Run fraction to boundary computation
a = run_fraction_to_boundary(p,q,z,d);

% Set rho
rho = z.rho;

% Try full steps for trial rhos
if rho < z.rho_last

  % Set rho_temp
  rho_temp = z.rho_last;

  % Loop
  while rho_temp > rho

    % Set rho
    z.rho = rho_temp;

    % Evaluate merit
    z = eval_merit(q,z);

    % Store current values
    phi = z.phi; x = z.x;
    if q.nE > 0, r1 = z.r1; r2 = z.r2; lambdaE = z.lambdaE; end;
    if q.nI > 0, s1 = z.s1; s2 = z.s2; lambdaI = z.lambdaI; end;

    % Try full step
    [c,z] = update_iterate(c,p,q,z,d,a,0);

    % Check for fraction-to-boundary violation
    ftb = 0;
    if q.nE > 0, ftb = ftb + sum(z.r1<p.ls_frac*r1) + sum(z.r2<p.ls_frac*r2); end;
    if q.nI > 0, ftb = ftb + sum(z.s1<p.ls_frac*s1) + sum(z.s2<p.ls_frac*s2); end;

    % Reset variables
    z.x = x;
    if q.nE > 0, z.r1 = r1; z.r2 = r2; z.lambdaE = lambdaE; end;
    if q.nI > 0, z.s1 = s1; z.s2 = s2; z.lambdaI = lambdaI; end;

    % Check Armijo condition
    if ftb == 0 & z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0), return; end;

    % Decrease rho
    rho_temp = p.rho_factor*rho_temp;
  
  end

end

% Set rho
z.rho = rho;

% Evaluate merit
z = eval_merit(q,z);

% Store current values
phi = z.phi; x = z.x;
if q.nE > 0, r1 = z.r1; r2 = z.r2; lambdaE = z.lambdaE; end;
if q.nI > 0, s1 = z.s1; s2 = z.s2; lambdaI = z.lambdaI; end;

% Backtracking loop
while a.p >= eps

  % Set trial point (this moves x,r1,r2,s1,r2,lambdaE,lambdaI)
  [c,z] = update_iterate(c,p,q,z,d,a,0);
  
  % Check for fraction-to-boundary violation
  ftb = 0;
  if q.nE > 0, ftb = ftb + sum(z.r1<p.ls_frac*r1) + sum(z.r2<p.ls_frac*r2); end;
  if q.nI > 0, ftb = ftb + sum(z.s1<p.ls_frac*s1) + sum(z.s2<p.ls_frac*s2); end;
  
  % Reset variables
  z.x = x;
  if q.nE > 0, z.r1 = r1; z.r2 = r2; z.lambdaE = lambdaE; end;
  if q.nI > 0, z.s1 = s1; z.s2 = s2; z.lambdaI = lambdaI; end;

  % Check Armijo condition
  if ftb == 0 & z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0), return; end;
  
  % Reduce steplength
  a.p = p.ls_factor*a.p;
  
end
