function x = eval_x_original(q,z)

% function x = eval_x_original(q,z)
%
% Author       : Frank E. Curtis
% Description  : Evaluates primal iterate in original space.
% Input        : q ~ quantities
%                z ~ iterate
% Output       : x ~ primal iterate in original space
% Last revised : 21 June 2010

% Construct x in original space
x       = zeros(q.nV_orig,1);
x(q.I1) = z.x(1               :q.n1               );
x(q.I2) = q.b2;
x(q.I3) = z.x(1+q.n1          :q.n1+q.n3          );
x(q.I4) = z.x(1+q.n1+q.n3     :q.n1+q.n3+q.n4     );
x(q.I5) = z.x(1+q.n1+q.n3+q.n4:q.n1+q.n3+q.n4+q.n5);
