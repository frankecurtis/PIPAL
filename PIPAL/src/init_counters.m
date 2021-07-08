function c = init_counters

% function c = init_counters
%
% Author       : Frank E. Curtis
% Description  : Initializes counters for iterations, function evaluations,
%                gradient evaluations, Hessian evaluations, and matrix
%                factorizations; initializes clock.
% Output       : c ~ counters
% Last revised : 21 June 2010

% Initialize iteration counter
c.k = 0;

% Initialize evaluation counters
c.f = 0;
c.g = 0;
c.H = 0;

% Initialize linear system solve counter
c.fact = 0;

% Start clock
tic;
