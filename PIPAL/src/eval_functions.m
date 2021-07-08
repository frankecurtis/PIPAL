function [c,z] = eval_functions(c,q,z)

% function [c,z] = eval_functions(c,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates objective and constraint values.
% Input        : c ~ counters
%                q ~ quantities
%                z ~ iterate
% Output       : c ~ updated counters
%                z ~ updated iterate
% Last revised : 21 June 2010

% Evaluate x in original space
x_orig = eval_x_original(q,z);

% Evaluate objective and constraints in AMPL
[z.f,c_orig] = spamfunc(x_orig,0);

% Increment function evaluation counter
c.f = c.f + 1;

% Set equality constraint values
if q.nE > 0, z.cE = c_orig(q.I6) - q.b6; end;

% Initialize inequality constraint values
if q.nI > 0, z.cI = zeros(q.nI,1); end;

% Set inequality constraint values
if q.n3 > 0, z.cI(1                                   :q.n3                                   ) =  q.l3 - z.x(1+q.n1          :q.n1+q.n3          ); end;
if q.n4 > 0, z.cI(1+q.n3                              :q.n3+q.n4                              ) = -q.u4 + z.x(1+q.n1+q.n3     :q.n1+q.n3+q.n4     ); end;
if q.n5 > 0, z.cI(1+q.n3+q.n4                         :q.n3+q.n4+q.n5                         ) =  q.l5 - z.x(1+q.n1+q.n3+q.n4:q.n1+q.n3+q.n4+q.n5);
             z.cI(1+q.n3+q.n4+q.n5                    :q.n3+q.n4+q.n5+q.n5                    ) = -q.u5 + z.x(1+q.n1+q.n3+q.n4:q.n1+q.n3+q.n4+q.n5); end;
if q.n7 > 0, z.cI(1+q.n3+q.n4+q.n5+q.n5               :q.n3+q.n4+q.n5+q.n5+q.n7               ) =  q.l7 - c_orig(q.I7);                              end;
if q.n8 > 0, z.cI(1+q.n3+q.n4+q.n5+q.n5+q.n7          :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8          ) = -q.u8 + c_orig(q.I8);                              end;
if q.n9 > 0, z.cI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8     :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9     ) =  q.l9 - c_orig(q.I9);
             z.cI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9:q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9+q.n9) = -q.u9 + c_orig(q.I9);                              end;

% Scale objective
z.f = z.f_scale*z.f;

% Scale equality constraints
if q.nE > 0, z.cE = z.cE_scale.*z.cE; end;

% Scale inequality constraints
if q.nI > 0, z.cI = z.cI_scale.*z.cI; end;
