function z = eval_slacks(p,q,z)

% function z = eval_slacks(p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates slack variables.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Check for equality constraints
if q.nE > 0

  % Set r1
  z.r1 = (1/2)*(z.mu - z.cE + sqrt(z.cE.^2 + z.mu^2));

  % Set r2
  z.r2 = (1/2)*(z.mu + z.cE + sqrt(z.cE.^2 + z.mu^2));

  % Adjust for numerical error
  z.r1 = max(z.r1,p.slack_min);
  z.r2 = max(z.r2,p.slack_min);

end

% Check for inequality constraints
if q.nI > 0

  % Set s1
  z.s1 = (1/2)*(2*z.mu - z.cI + sqrt(z.cI.^2 + 4*z.mu^2));

  % Set s2
  z.s2 = (1/2)*(2*z.mu + z.cI + sqrt(z.cI.^2 + 4*z.mu^2));

  % Adjust for numerical error
  z.s1 = max(z.s1,p.slack_min);
  z.s2 = max(z.s2,p.slack_min);

end
