function z = update_point(q,z,d,a)

% function z = update_point(q,z,d,a)
%
% Author       : Frank E. Curtis
% Description  : Updates primal and dual points.
% Input        : q ~ quantities
%                z ~ iterate
%                d ~ direction
%                a ~ step acceptance values
% Output       : z ~ updated iterate
% Last revised : 21 June 2010

% Cut length of step
d = run_cut_step(q,d,a);

% Update primal and dual variables
             z.x       = z.x       + d.x;
if q.nE > 0, z.r1      = z.r1      + d.r1;
             z.r2      = z.r2      + d.r2;      end;
if q.nI > 0, z.s1      = z.s1      + d.s1;
             z.s2      = z.s2      + d.s2;      end;
if q.nE > 0, z.lambdaE = z.lambdaE + d.lambdaE; end;
if q.nI > 0, z.lambdaI = z.lambdaI + d.lambdaI; end;
