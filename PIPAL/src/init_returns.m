function r = init_returns(p,q,z)

% function r = init_returns(p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Initializes returns including objective value,
%                feasibility violation, optimality error, and
%                primal and dual iterates.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : r ~ returns
% Last revised : 21 June 2010

% Initialize termination message
r.msg = '---';

% Initialize objective value
r.f = z.f;

% Initialize feasibility violation
r.v = z.v;

% Initialize optimality error
r.kkt = z.kkt;

% Initialize primal iterate
r.x = z.x;

% Initialize dual iterate for equality constraints
if q.nE > 0, r.lambdaE = z.lambdaE; end;

% Initialize dual iterate for inequality constraints
if q.nI > 0, r.lambdaI = z.lambdaI; end;
