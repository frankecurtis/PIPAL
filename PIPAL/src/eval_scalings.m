function z = set_scalings(p,q,z)

% function z = set_scalings(p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates objective and constraint scaling factors
%                to limit gradient norms at initial point.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : z ~ udpated iterate
% Last revised : 21 June 2010

% Evaluate norm of objective gradient
g_norm_inf = norm(z.g,inf);
  
% Scale down objective if norm of gradient is too large
z.f_scale = p.grad_max/max(g_norm_inf,p.grad_max);

% Initialize equality constraint scalings
if q.nE > 0, z.cE_scale = zeros(q.nE,1); end;

% Loop through equality constraints
for i = 1:q.nE
      
  % Evaluate norm of gradient of ith equality constraint
  JE_i_norm_inf = norm(z.JE(i,:),inf);

  % Scale down equality constraint i if norm of gradient is too large
  z.cE_scale(i) = p.grad_max/max(JE_i_norm_inf,p.grad_max);

end

% Initialize inequality constraint scalings
if q.nI > 0, z.cI_scale = zeros(q.nI,1); end;

% Loop through inequality constraints
for i = 1:q.nI
      
  % Evaluate norm of gradient of ith inequality constraint
  JI_i_norm_inf = norm(z.JI(i,:),inf);

  % Scale down inequality constraint i if norm of gradient is too large
  z.cI_scale(i) = p.grad_max/max(JI_i_norm_inf,p.grad_max);

end
