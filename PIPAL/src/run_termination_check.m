function b = run_termination_check(r)

% function b = run_termination_check(r)
%
% Author       : Frank E. Curtis
% Description  : Checks termination of iteration loop.
% Input        : r ~ returns
% Output       : b ~ 0 if done; 1 otherwise
% Last revised : 21 June 2010

% Check for termination
if strcmp(r.msg,'---') == 1, b = 1; else b = 0; end;
