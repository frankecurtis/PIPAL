% README.m

% Please cite:
%   Frank E. Curtis.  "A Penalty-Interior-Point Algorithm for Nonlinear Constrained
%   Optimization."  Mathematical Programming Computation, 4(2):181–209, 2012.

% This code runs a penalty-interior-point algorithm for solving nonlinear
% constrained optimization problems of the form
%   minimize f(x) subject to cE(x) = 0, cI(x) <= 0,
% where f, cE, and cI are assumed to be continuously differentiable on R^n.

% Usage: Let nl be the path of an AMPL .nl file and let algorithm be chosen
% in {0,1}.  Then, the AMPL model can be solved with the commands:
%   >> P = Pipal(nl,algorithm);
%   >> P.optimize;
% (An amplfunc.mex file that reads AMPL .nl files is necessary; please see
% D. Gay, ``Hooking Your Solver to AMPL,'' Technical Report, Computing
% Sciences Research Center, Bell Laboratories, Murray Hill, NJ, USA.)
% Algorithm option 0 runs a conservative updating strategy for the penalty
% and interior-point parameters and algorithm option 1 runs an aggressive
% strategy.  An output file pipal.out will be generated.