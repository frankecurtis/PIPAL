function z = eval_infeasibility(q,z)

% function z = eval_infeasibility(q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates feasibility violation.
% Input        : q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Initialize feasibility violation vector
vec = [];

% Check for equality constraints
if q.nE > 0, vec = z.cE; end;

% Check for inequality constraints
if q.nI > 0, vec = [vec; max(z.cI,0)]; end;

% Evaluate feasibility violation
z.v = norm(vec,1);
