function [c,z] = eval_gradients(c,q,z,opt)

% function [c,z] = eval_gradients(c,q,z,opt)
%
% Author       : Frank E. Curtis
% Description  : Evaluates objective and constraint gradients.
% Input        : c   ~ counters
%                q   ~ quantities
%                z   ~ iterate
%                opt ~ evaluation option
%                        0 ~ only computing
%                        1 ~ only scaling
%                        2 ~ computing and scaling
% Output       : c   ~ update counters
%                z   ~ updated iterate
% Last revised : 21 June 2010

% Check evaluation option
if opt ~= 1

  % Evaluate x in original space
  x_orig = eval_x_original(q,z);

  % Evaluate gradient and Jacobian in AMPL
  [g_orig,J_orig] = spamfunc(x_orig,1);

  % Increment gradient evaluation counter
  c.g = c.g + 1;

  % Set objective gradient
  z.g = [g_orig(q.I1); g_orig(q.I3); g_orig(q.I4); g_orig(q.I5)];
  
  % Set equality constraint Jacobian
  if q.nE > 0, z.JE = [J_orig(q.I6,q.I1) J_orig(q.I6,q.I3) J_orig(q.I6,q.I4) J_orig(q.I6,q.I5)]; end;

  % Initialize inequality constraint Jacobian
  if q.nI > 0, z.JI = sparse(q.nI,q.nV); end;

  % Set inequality constraint Jacobian
  if q.n3 > 0, z.JI(1                                   :q.n3                                   ,1+q.n1          :q.n1+q.n3          ) = -speye(q.n3); end;
  if q.n4 > 0, z.JI(1+q.n3                              :q.n3+q.n4                              ,1+q.n1+q.n3     :q.n1+q.n3+q.n4     ) =  speye(q.n4); end;
  if q.n5 > 0, z.JI(1+q.n3+q.n4                         :q.n3+q.n4+q.n5                         ,1+q.n1+q.n3+q.n4:q.n1+q.n3+q.n4+q.n5) = -speye(q.n5);
               z.JI(1+q.n3+q.n4+q.n5                    :q.n3+q.n4+q.n5+q.n5                    ,1+q.n1+q.n3+q.n4:q.n1+q.n3+q.n4+q.n5) =  speye(q.n5); end;
  if q.n7 > 0, z.JI(1+q.n3+q.n4+q.n5+q.n5               :q.n3+q.n4+q.n5+q.n5+q.n7               ,1               :q.n1+q.n3+q.n4+q.n5) = -[J_orig(q.I7,q.I1) J_orig(q.I7,q.I3) J_orig(q.I7,q.I4) J_orig(q.I7,q.I5)]; end;
  if q.n8 > 0, z.JI(1+q.n3+q.n4+q.n5+q.n5+q.n7          :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8          ,1               :q.n1+q.n3+q.n4+q.n5) =  [J_orig(q.I8,q.I1) J_orig(q.I8,q.I3) J_orig(q.I8,q.I4) J_orig(q.I8,q.I5)]; end;
  if q.n9 > 0, z.JI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8     :q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9     ,1               :q.n1+q.n3+q.n4+q.n5) = -[J_orig(q.I9,q.I1) J_orig(q.I9,q.I3) J_orig(q.I9,q.I4) J_orig(q.I9,q.I5)];
               z.JI(1+q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9:q.n3+q.n4+q.n5+q.n5+q.n7+q.n8+q.n9+q.n9,1               :q.n1+q.n3+q.n4+q.n5) =  [J_orig(q.I9,q.I1) J_orig(q.I9,q.I3) J_orig(q.I9,q.I4) J_orig(q.I9,q.I5)]; end;

end

% Check evaluation option
if opt > 0

  % Scale objective gradient
  z.g = z.f_scale*z.g;

  % Scale equality constraint Jacobian
  if q.nE > 0, z.JE = spdiags(z.cE_scale,0:0,q.nE,q.nE)*z.JE; end;
  
  % Scale inequality constraint Jacobian
  if q.nI > 0, z.JI = spdiags(z.cI_scale,0:0,q.nI,q.nI)*z.JI; end;

end
