function z = eval_merit(q,z)

% function z = eval_merit(q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates penalty-interior-point objective; i.e., merit function.
% Input        : q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Initialize penalty-interior-point objective
z.phi = z.rho*z.f;

% Check for equality constraints
if q.nE > 0, z.phi = z.phi - z.mu*sum(log([z.r1;z.r2])) + sum([z.r1;z.r2]); end;

% Check for inequality constraints
if q.nI > 0, z.phi = z.phi - z.mu*sum(log([z.s1;z.s2])) + sum(z.s2);        end;
