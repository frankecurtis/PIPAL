function print_step_acceptance(o,a)

% function print_step_acceptance(o,a)
%
% Author       : Frank E. Curtis
% Description  : Prints step acceptance information.
% Input        : o ~ outputs
%                a ~ step acceptances
% Last revised : 21 June 2010

% Print iterate information
if o.verbosity <= 0
  fprintf(o.fout,'%.4e\n',a.p);
else
  fprintf(o.fout,'%.4e  %.4e\n',a.p,a.d);
end
