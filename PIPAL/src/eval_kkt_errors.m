function z = eval_kkt_errors(q,z)

% function z = eval_kkt_errors(q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates optimality errors:
%                  z.kkt(1) ~ corresponds to rho =   0, mu =  0
%                  z.kkt(2) ~ corresponds to rho = rho, mu =  0
%                  z.kkt(3) ~ corresponds to rho = rho, mu = mu
% Input        : q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Loop to compute optimality errors
z.kkt(1) = eval_kkt_error(q,z,0    ,0   );
z.kkt(2) = eval_kkt_error(q,z,z.rho,0   );
z.kkt(3) = eval_kkt_error(q,z,z.rho,z.mu);

% Error function
function kkt_error = eval_kkt_error(q,z,rho,mu)

% Initialize optimality vector
kkt = zeros(q.nV+2*q.nE+2*q.nI,1);

% Set gradient of penalty objective
kkt(1:q.nV) = rho*z.g;

% Set gradient of Lagrangian for equality constraints
if q.nE > 0, kkt(1:q.nV) = kkt(1:q.nV) + (z.lambdaE'*z.JE)'; end;

% Set gradient of Lagrangian for inequality constraints
if q.nI > 0, kkt(1:q.nV) = kkt(1:q.nV) + (z.lambdaI'*z.JI)'; end;

% Set complementarity for equality constraint slacks
if q.nE > 0, kkt(1+q.nV:q.nV+2*q.nE) = [z.r1.*(1 + z.lambdaE) - mu;
                                        z.r2.*(1 - z.lambdaE) - mu]; end;

% Set complementarity for inequality constraint slacks
if q.nI > 0, kkt(1+q.nV+2*q.nE:q.nV+2*q.nE+2*q.nI) = [z.s1.*(0 + z.lambdaI) - mu;
                                                      z.s2.*(1 - z.lambdaI) - mu]; end;

% Scale complementarity
if rho > 0, kkt = (1/max(1,norm(rho*z.g,inf)))*kkt; end;

% Evaluate optimality error
kkt_error = norm(kkt,inf);
