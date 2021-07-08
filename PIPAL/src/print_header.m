function print_header(o,q,z)

% function print_header(o,q,z)
%
% Author       : Frank E. Curtis
% Description  : Prints output header.
% Input        : o ~ outputs
%                q ~ quantities
%                z ~ iterate
% Last revised : 21 June 2010

% Print problem size
fprintf(o.fout,'Problem size\n');
fprintf(o.fout,'============\n');
fprintf(o.fout,'  Number of variables....................... : %8d\n',q.nV);
fprintf(o.fout,'  Number of equality constraints............ : %8d\n',q.nE);
fprintf(o.fout,'  Number of inequality constraints.......... : %8d\n',q.nI);
fprintf(o.fout,'\n');

% Print problem nonzeros
fprintf(o.fout,'Problem sparsity\n');
fprintf(o.fout,'================\n');
fprintf(o.fout,'  Nonzeros in Hessian of Lagrangian......... : %8d\n',nnz(z.H));
if q.nE > 0, fprintf(o.fout,'  Nonzeros in equality constraint Jacobian.. : %8d\n',nnz(z.JE));
else         fprintf(o.fout,'  Nonzeros in equality constraint Jacobian.. : %8d\n',0); end;
if q.nI > 0, fprintf(o.fout,'  Nonzeros in inequality constraint Jacobian : %8d\n',nnz(z.JI));
else         fprintf(o.fout,'  Nonzeros in inequality constraint Jacobian : %8d\n',0); end;
fprintf(o.fout,'\n');
