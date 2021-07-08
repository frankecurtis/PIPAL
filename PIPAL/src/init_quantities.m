function q = init_quantities(i)

% function q = init_quantities(i)
%
% Author       : Frank E. Curtis
% Description  : Initializes AMPL-derived quantities; determines
%                equality and inequality constraints; determines
%                problem sizes; determines initial point.
% Input        : i ~ inputs
% Output       : q ~ quantities
% Last revised : 21 June 2010

% Set initial AMPL data
[x_orig,x_lo,x_hi,l_orig,c_lo,c_hi] = spamfunc(i.nl);

% Find index sets
q.I1 = find(x_lo == -inf & x_hi == inf);
q.I2 = find(x_lo == x_hi);
q.I3 = find(x_lo >  -inf & x_hi == inf);
q.I4 = find(x_lo == -inf & x_hi <  inf);
q.I5 = find(x_lo >  -inf & x_hi <  inf & x_lo ~= x_hi);
q.I6 = find(c_lo == c_hi);
q.I7 = find(c_lo >  -inf & c_hi == inf);
q.I8 = find(c_lo == -inf & c_hi <  inf);
q.I9 = find(c_lo >  -inf & c_hi <  inf & c_lo ~= c_hi);

% Set rhs values
q.b2 = x_lo(q.I2);
q.l3 = x_lo(q.I3);
q.u4 = x_hi(q.I4);
q.l5 = x_lo(q.I5);
q.u5 = x_hi(q.I5);
q.b6 = c_lo(q.I6);
q.l7 = c_lo(q.I7);
q.u8 = c_hi(q.I8);
q.l9 = c_lo(q.I9);
q.u9 = c_hi(q.I9);

% Set sizes of index sets
q.n1 = length(q.I1);
q.n2 = length(q.I2);
q.n3 = length(q.I3);
q.n4 = length(q.I4);
q.n5 = length(q.I5);
q.n6 = length(q.I6);
q.n7 = length(q.I7);
q.n8 = length(q.I8);
q.n9 = length(q.I9);

% Set number of variables and constraints
q.nV = q.n1 + q.n3 + q.n4 + q.n5;
q.nI = q.n3 + q.n4 + 2*q.n5 + q.n7 + q.n8 + 2*q.n9;
q.nE = q.n6;

% Set number of original variables
q.nV_orig = q.n1 + q.n2 + q.n3 + q.n4 + q.n5;

% Set size of primal-dual matrix
q.size = q.nV + 3*q.nE + 3*q.nI;

% Set initial point
q.x_0 = [x_orig(q.I1); x_orig(q.I3); x_orig(q.I4); x_orig(q.I5);];
