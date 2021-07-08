function d = run_kkt_solver(p,q,z)

% function d = run_kkt_solver(p,q,z)
%
% Author       : Frank E. Curtis
% Description  : Runs primal-dual system solver.
% Input        : p ~ parameters
%                q ~ quantities
%                z ~ iterate
% Output       : d ~ direction
% Last revised : 21 June 2010

% Evaluate direction
direction = z.kkt_matrix_S(:,z.kkt_matrix_P)*(z.kkt_matrix_L'\(z.kkt_matrix_D\(z.kkt_matrix_L\(z.kkt_matrix_S(z.kkt_matrix_P,:)*(-z.kkt_rhs)))));

% Parse direction
             d.x       = direction(1                              :q.nV                              );
if q.nE > 0, d.r1      = direction(1+q.nV                         :q.nV+q.nE                         );
             d.r2      = direction(1+q.nV+q.nE                    :q.nV+q.nE+q.nE                    ); end;
if q.nI > 0, d.s1      = direction(1+q.nV+q.nE+q.nE               :q.nV+q.nE+q.nE+q.nI               );
             d.s2      = direction(1+q.nV+q.nE+q.nE+q.nI          :q.nV+q.nE+q.nE+q.nI+q.nI          ); end;
if q.nE > 0, d.lambdaE = direction(1+q.nV+q.nE+q.nE+q.nI+q.nI     :q.nV+q.nE+q.nE+q.nI+q.nI+q.nE     ); end;
if q.nI > 0, d.lambdaI = direction(1+q.nV+q.nE+q.nE+q.nI+q.nI+q.nE:q.nV+q.nE+q.nE+q.nI+q.nI+q.nE+q.nI); end;

% Evaluate primal direction norm
d.x_norm = norm(d.x);

% Initialize dual direction
lambda = [];

% Update dual direction
if q.nE > 0, lambda = d.lambdaE; end;
if q.nI > 0, lambda = [lambda; d.lambdaI]; end;

% Evaluate dual direction norm
d.lambda_norm = norm(lambda);
