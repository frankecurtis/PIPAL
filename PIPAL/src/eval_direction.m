function d = eval_direction(q,d1,d2,d3,a1,a2,a3,z)

% function d = eval_direction(q,d1,d2,d3,a1,a2,a3,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates direction from linear combination of d1, d2, and d3.
% Input        : q  ~ quantities
%                d1 ~ direction
%                d2 ~ direction
%                d3 ~ direction
%                a1 ~ scalar
%                a2 ~ scalar
%                a3 ~ scalar
%                z  ~ iterate
% Output       : d  ~ a1*d1 + a2*d2 + a3*d3
% Last revised : 21 June 2010

% Evaluate linear combinations
             d.x       = a1*d1.x       + a2*d2.x       + a3*d3.x;
if q.nE > 0, d.r1      = a1*d1.r1      + a2*d2.r1      + a3*d3.r1;
             d.r2      = a1*d1.r2      + a2*d2.r2      + a3*d3.r2;      end;
if q.nI > 0, d.s1      = a1*d1.s1      + a2*d2.s1      + a3*d3.s1;
             d.s2      = a1*d1.s2      + a2*d2.s2      + a3*d3.s2;      end;
if q.nE > 0, d.lambdaE = a1*d1.lambdaE + a2*d2.lambdaE + a3*d3.lambdaE; end;
if q.nI > 0, d.lambdaI = a1*d1.lambdaI + a2*d2.lambdaI + a3*d3.lambdaI; end;

% Set primal direction norm
d.x_norm = norm(d.x);

% Initialize dual direction
lambda = [];

% Update dual direction
if q.nE > 0, lambda = d.lambdaE; end;
if q.nI > 0, lambda = [lambda; d.lambdaI]; end;

% Evaluate dual direction norm
d.lambda_norm = norm(lambda);
