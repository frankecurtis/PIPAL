function print_direction(o,z,d)

% function print_direction(o,z,d)
%
% Author       : Frank E. Curtis
% Description  : Prints direction information.
% Input        : o ~ outputs
%                z ~ iterate
%                d ~ direction
% Last revised : 21 June 2010

% Print step information
if o.verbosity <= 0
  fprintf(o.fout,'%.4e | ',d.x_norm);
else
  fprintf(o.fout,'%+.4e  %.4e | %.4e  %.4e  %.4e  %+.4e  %+.4e  %+.4e | ',z.phi,z.kkt(3),z.shift,d.x_norm,d.lambda_norm,d.ltred0,d.qtred,d.m);
end
