function d = eval_models(q,z,d)

% function d = eval_models(q,z,d)
%
% Author       : Frank E. Curtis
% Description  : Evaluates models of merit function along direction.
% Input        : q ~ quantities
%                z ~ iterate
%                d ~ direction
% Output       : d ~ updated direction
% Last revised : 21 June 2010

% Evaluate linear model of penalty-interior-point objective for zero penalty parameter
d.lred0 = 0;
if q.nE > 0, d.lred0 = d.lred0 - sum([1-z.mu./z.r1; 1-z.mu./z.r2].*[d.r1; d.r2]); end;
if q.nI > 0, d.lred0 = d.lred0 - sum([0-z.mu./z.s1; 1-z.mu./z.s2].*[d.s1; d.s2]); end;

% Evaluate remaining quantities only for nonzero penalty parameter
if z.rho > 0

  % Evaluate linear model of merit function for zero penalty parameter
  d.ltred0 = 0;
  if q.nE > 0, d.ltred0 = d.ltred0 - (1/2)*full(sum(((1-z.mu./z.r1).*(-1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))+(1-z.mu./z.r2).*(1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))).*(z.JE*d.x))); end;
  if q.nI > 0, d.ltred0 = d.ltred0 - (1/2)*full(sum(((0-z.mu./z.s1).*(-1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))+(1-z.mu./z.s2).*(1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))).*(z.JI*d.x))); end;

  % Evaluate linear model of merit function
  d.ltred = -z.rho*z.g'*d.x + d.ltred0;

  % Evaluate quadratic model of merit function
  d.qtred = d.ltred - (1/2)*d.x'*z.H*d.x;
  if q.nE > 0, Jd = z.JE*d.x; Dinv = z.r1./(1+z.lambdaE)+z.r2./(1-z.lambdaE); d.qtred = d.qtred - (1/2)*Jd'*(Jd./Dinv); end;
  if q.nI > 0, Jd = z.JI*d.x; Dinv = z.s1./(0+z.lambdaI)+z.s2./(1-z.lambdaI); d.qtred = d.qtred - (1/2)*Jd'*(Jd./Dinv); end;

  % Initialize quality function vector
  vec = zeros(q.nV+2*q.nE+2*q.nI,1);

  % Set gradient of objective
  vec(1:q.nV) = z.rho*z.g;

  % Set gradient of Lagrangian for equality constraints
  if q.nE > 0, vec(1:q.nV) = vec(1:q.nV) + ((z.lambdaE+d.lambdaE)'*z.JE)'; end;

  % Set gradient of Lagrangian for inequality constraints
  if q.nI > 0, vec(1:q.nV) = vec(1:q.nV) + ((z.lambdaI+d.lambdaI)'*z.JI)'; end;

  % Set complementarity for equality constraint slacks
  if q.nE > 0, vec(1+q.nV:q.nV+2*q.nE) = [(z.r1+d.r1).*(1 + (z.lambdaE+d.lambdaE));
                                          (z.r2+d.r2).*(1 - (z.lambdaE+d.lambdaE))]; end;

  % Set complementarity for inequality constraint slacks
  if q.nI > 0, vec(1+q.nV+2*q.nE:q.nV+2*q.nE+2*q.nI) = [(z.s1+d.s1).*(0 + (z.lambdaI+d.lambdaI));
                                                        (z.s2+d.s2).*(1 - (z.lambdaI+d.lambdaI))]; end;

  % Evaluate quality function
  d.m = norm(vec,inf);

end
