function r = update_returns(c,p,q,z,r)

% function r = update_returns(c,p,q,z,r)
%
% Author       : Frank E. Curtis
% Description  : Updates returns; checks termination tolerances.
% Input        : c ~ counters
%                p ~ parameters
%                q ~ quantities
%                z ~ iterate
%                r ~ returns
% Output       : r ~ updated returns
% Last revised : 21 June 2010

% Update termination based on iteration count
if c.k >= p.iter_max, r.msg = 'itr'; end;

% Update termination based on optimality error of nonlinear optimization problem
if z.kkt(2) <= p.opt_err_tol & z.v <= p.opt_err_tol*max(1,z.v_0), r.msg = 'opt'; end;

% Update termination based on optimality error of feasibility problem
if z.kkt(1) <= p.opt_err_tol & z.v >  p.opt_err_tol*max(1,z.v_0), r.msg = 'inf'; end;

% Update objective
r.f = z.f;

% Update feasibility violation
r.v = z.v;

% Update optimality error
r.kkt = z.kkt;

% Update primal iterate
r.x = z.x;

% Update dual iterate for equality constraints
if q.nE > 0, r.lambdaE = z.lambdaE; end;

% Update dual iterate for inequality constraints
if q.nI > 0, r.lambdaI = z.lambdaI; end;
