function r = run_driver(i)

% function r = run_driver(i)
%
% Author       : Frank E. Curtis
% Description  : Runs driver for optimization; verifies
%                AMPL file is provided and exists.
% Input        : i ~ inputs
% Output       : r ~ returns
% Last revised : 21 June 2010

% Assert that AMPL file has been provided
assert(isfield(i,'nl'),'PIPAL: AMPL file, i.nl, not specified.');

% Assert that AMPL file exists
assert(exist(i.nl,'file')~=0,sprintf('PIPAL: AMPL file, %s, does not exist.\n',i.nl));

% Run optimization algorithm
r = run_optimization(i);
