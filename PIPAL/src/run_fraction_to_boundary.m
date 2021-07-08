function a = run_fraction_to_boundary(p,q,z,d)

% function a = run_fraction_to_boundary(p,q,z,d)
%
% Author       : Frank E. Curtis
% Description  : Runs fraction to boundary computation.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
%                d ~ direction
% Output       : a ~ step acceptance values
% Last revised : 21 June 2010

% Initialize primal fraction-to-boundary
a.p0 = 1;

% Update primal fraction-to-boundary for equality constraint slacks
if q.nE > 0, a.p0 = min([a.p0; (p.ls_frac-1)*z.r1(d.r1<0)./d.r1(d.r1<0); (p.ls_frac-1)*z.r2(d.r2<0)./d.r2(d.r2<0)]); end;

% Update primal fraction-to-boundary for inequality constraint slacks
if q.nI > 0, a.p0 = min([a.p0; (p.ls_frac-1)*z.s1(d.s1<0)./d.s1(d.s1<0); (p.ls_frac-1)*z.s2(d.s2<0)./d.s2(d.s2<0)]); end;

% Initialize primal steplength
a.p = a.p0;

% Initialize dual fraction-to-boundary
a.d = 1;

% Update dual fraction-to-boundary for equality constraint multipliers
if q.nE > 0, a.d = min([a.d; (p.ls_frac-1)*(1+z.lambdaE(d.lambdaE<0))./d.lambdaE(d.lambdaE<0); (1-p.ls_frac)*(1-z.lambdaE(d.lambdaE>0))./d.lambdaE(d.lambdaE>0)]); end;

% Update dual fraction-to-boundary for inequality constraint multipliers
if q.nI > 0, a.d = min([a.d; (p.ls_frac-1)*(0+z.lambdaI(d.lambdaI<0))./d.lambdaI(d.lambdaI<0); (1-p.ls_frac)*(1-z.lambdaI(d.lambdaI>0))./d.lambdaI(d.lambdaI>0)]); end;
