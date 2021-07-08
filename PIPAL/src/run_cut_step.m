function d = run_cut_step(q,d,a)

% function d = run_cut_step(q,d,a)
%
% Author       : Frank E. Curtis
% Description  : Runs routine to cut length of step.
% Input        : q ~ quantities
%                d ~ direction
%                a ~ step acceptance values
% Output       : d ~ updated direction
% Last revised : 21 June 2010

% Cut direction
             d.x       = a.p*d.x;
if q.nE > 0, d.r1      = a.p*d.r1;
             d.r2      = a.p*d.r2;      end;
if q.nI > 0, d.s1      = a.p*d.s1;
             d.s2      = a.p*d.s2;      end;
if q.nE > 0, d.lambdaE = a.d*d.lambdaE; end;
if q.nI > 0, d.lambdaI = a.d*d.lambdaI; end;
