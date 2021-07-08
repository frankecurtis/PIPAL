function [c,z] = eval_hessian(c,q,z)

% function [c,z] = eval_hessian(c,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates Hessian of Lagrangian.
% Input        : c ~ counters
%                q ~ quantities
%                z ~ iterate
% Output       : c ~ updated counters
%                z ~ updated iterate
% Last revised : 21 June 2010

% Initialize multipliers in original space
l_orig = zeros(q.nE+q.n7+q.n8+q.n9,1);

% Scale equality constraint multipliers
if q.nE > 0, lambdaE = z.lambdaE.*(z.cE_scale/(z.rho*z.f_scale)); end;

% Set equality constraint multipliers in original space
if q.nE > 0, l_orig(q.I6) = lambdaE; end;

% Scale inequality constraint multipliers
if q.n7+q.n8+q.n9 > 0, lambdaI = z.lambdaI.*(z.cI_scale/(z.rho*z.f_scale)); end;

% Set inequality constraint multipliers in original space
if q.n7 > 0, l_orig(q.I7) = -lambdaI(1+q.n3+q.n4+q.n5+q.n5               :q.n3+q.n4+q.n5+q.n5+q.n7               ); end;
if q.n8 > 0, l_orig(q.I8) = +lambdaI(1+q.n3+q.n4+q.n5+q.n5+q.n7          :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8          ); end;
if q.n9 > 0, l_orig(q.I9) = -lambdaI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8     :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9     )...
                            +lambdaI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9:q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9+q.n9); end;

% Evaluate H_orig
if (q.nE+q.n7+q.n8+q.n9 == 0), H_orig = spamfunc([]);
else,                          H_orig = spamfunc(l_orig); end;

% Increment Hessian evaluation counter
c.H = c.H + 1;

% Set Hessian of the Lagrangian
z.H = [H_orig(q.I1,q.I1) H_orig(q.I1,q.I3) H_orig(q.I1,q.I4) H_orig(q.I1,q.I5);
       H_orig(q.I3,q.I1) H_orig(q.I3,q.I3) H_orig(q.I3,q.I4) H_orig(q.I3,q.I5);
       H_orig(q.I4,q.I1) H_orig(q.I4,q.I3) H_orig(q.I4,q.I4) H_orig(q.I4,q.I5);
       H_orig(q.I5,q.I1) H_orig(q.I5,q.I3) H_orig(q.I5,q.I4) H_orig(q.I5,q.I5);];

% Rescale H and Ht
z.H = z.rho*z.f_scale*z.H;
