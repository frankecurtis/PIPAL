function z = eval_kkt_rhs(q,z)

% function z = eval_kkt_rhs(q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates primal-dual system right-hand side.
% Input        : q ~ quantities
%                z ~ iterate
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Initialize optimality vector
z.kkt_rhs = zeros(q.size,1);

% Set gradient of objective
z.kkt_rhs(1:q.nV) = z.rho*z.g;

% Set gradient of Lagrangian for equality constraints
if q.nE > 0, z.kkt_rhs(1:q.nV) = z.kkt_rhs(1:q.nV) + (z.lambdaE'*z.JE)'; end;

% Set gradient of Lagrangian for inequality constraints
if q.nI > 0, z.kkt_rhs(1:q.nV) = z.kkt_rhs(1:q.nV) + (z.lambdaI'*z.JI)'; end;

% Set complementarity for equality constraint slacks
if q.nE > 0, z.kkt_rhs(1+q.nV:q.nV+2*q.nE) = [1 + z.lambdaE - z.mu./z.r1;
                                              1 - z.lambdaE - z.mu./z.r2]; end;

% Set complementarity for inequality constraint slacks
if q.nI > 0, z.kkt_rhs(1+q.nV+2*q.nE:q.nV+2*q.nE+2*q.nI) = [0 + z.lambdaI - z.mu./z.s1;
                                                            1 - z.lambdaI - z.mu./z.s2]; end;

% Set penalty-interior-point constraint values for equality constraints
if q.nE > 0, z.kkt_rhs(1+q.nV+2*q.nE+2*q.nI:q.nV+3*q.nE+2*q.nI) = z.cE + z.r1 - z.r2; end;

% Set penalty-interior-point constraint values for inequality constraints
if q.nI > 0, z.kkt_rhs(1+q.nV+3*q.nE+2*q.nI:q.nV+3*q.nE+3*q.nI) = z.cI + z.s1 - z.s2; end;
