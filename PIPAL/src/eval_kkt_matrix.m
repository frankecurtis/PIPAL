function [c,z] = eval_kkt_matrix(c,p,q,z)

% function [c,z] = eval_kkt_matrix(c,p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates primal-dual system matrix;
%                corrects inertia, if necessary.
% Input        : c ~ counters
%                p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : c ~ updated counters
%                z ~ updated iterate
% Last revised : 21 June 2010

% Set interior-point Hessian for equality constraints
if q.nE > 0, z.O = spdiags([(1+z.lambdaE)./z.r1; (1-z.lambdaE)./z.r2],0:0,2*q.nE,2*q.nE); end;

% Set interior-point Hessian for inequality constraints
if q.nI > 0, z.G = spdiags([(0+z.lambdaI)./z.s1; (1-z.lambdaI)./z.s2],0:0,2*q.nI,2*q.nI); end;

% Set minimum potential shift
if isfield(z,'shift'), min_shift = max(p.shift_min,p.shift_factor*z.shift); else min_shift = p.shift_min; end;

% Initialize Hessian modification
if isfield(z,'cut_last') & z.cut_last == 1, z.shift = min(p.shift_max,min_shift/p.shift_factor); else z.shift = 0; end;

% Initialize inertia correction loop
done = 0; z.shift22 = 0;

% Loop until inertia is correct
while ~done & z.shift < p.shift_max

  % Initialize primal-dual matrix
  z.kkt_matrix = sparse(q.size,q.size);
  
  % Set Hessian of Lagrangian
  z.kkt_matrix(1:q.nV,1:q.nV) = z.H+z.shift*speye(q.nV);
  
  % Check for equality constraints
  if q.nE > 0
  
    % Set barrier Hessian
    z.kkt_matrix(1+q.nV:q.nV+2*q.nE,1+q.nV:q.nV+2*q.nE) = z.O;
    
    % Set constraint Jacobian
    z.kkt_matrix(1+q.nV+2*q.nE+2*q.nI:q.nV+3*q.nE+2*q.nI,:) = [z.JE speye(q.nE) -speye(q.nE) sparse(q.nE,q.nI) sparse(q.nE,q.nI) -z.shift22*speye(q.nE) sparse(q.nE,q.nI)];
    z.kkt_matrix(:,1+q.nV+2*q.nE+2*q.nI:q.nV+3*q.nE+2*q.nI) = [z.JE speye(q.nE) -speye(q.nE) sparse(q.nE,q.nI) sparse(q.nE,q.nI) -z.shift22*speye(q.nE) sparse(q.nE,q.nI)]';
  
  end
  
  % Check for inequality constraints
  if q.nI > 0
  
    % Set barrier Hessian
    z.kkt_matrix(1+q.nV+2*q.nE:q.nV+2*q.nE+2*q.nI,1+q.nV+2*q.nE:q.nV+2*q.nE+2*q.nI) = z.G;

    % Set constraint Jacobian
    z.kkt_matrix(1+q.nV+3*q.nE+2*q.nI:q.nV+3*q.nE+3*q.nI,:) = [z.JI sparse(q.nI,q.nE) sparse(q.nI,q.nE) speye(q.nI) -speye(q.nI) sparse(q.nI,q.nE) -z.shift22*speye(q.nI)];    
    z.kkt_matrix(:,1+q.nV+3*q.nE+2*q.nI:q.nV+3*q.nE+3*q.nI) = [z.JI sparse(q.nI,q.nE) sparse(q.nI,q.nE) speye(q.nI) -speye(q.nI) sparse(q.nI,q.nE) -z.shift22*speye(q.nI)]';    
  
  end

  % Factor primal-dual matrix  
  [z.kkt_matrix_L,z.kkt_matrix_D,z.kkt_matrix_P,z.kkt_matrix_S,neig] = ldl(tril(z.kkt_matrix),p.pivot_thresh,'vector');
  
  % Increment factorization counter
  c.fact = c.fact + 1;

  % Set number of nonnegative eigenvalues
  peig = q.size - neig;
  
  % Check inertia
  if     peig < q.nV+2*q.nE+2*q.nI        , z.shift   = max(min_shift,z.shift/p.shift_factor);
  elseif neig < q.nE+q.nI & z.shift22 == 0, z.shift22 = p.shift_min;
  else                                    , done      = 1; end;

end

% Update Hessian
z.H = z.H+z.shift*speye(q.nV);
