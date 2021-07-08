function print_footer(o,c,z,r)

% function print_footer(o,c,z,r)
%
% Author       : Frank E. Curtis
% Description  : Prints output footer.
% Input        : o ~ outputs
%                c ~ counters
%                z ~ iterate
%                r ~ returns
% Last revised : 21 June 2010

% Print final iterate
print_iterate(o,c,z);

% Print close of algorithm output
fprintf(o.fout,'%s\n%s\n',o.none,o.line);
fprintf(o.fout,'\n');

% Print solver result
fprintf(o.fout,'Final result\n');
fprintf(o.fout,'============\n');
if strcmp(r.msg,'---') == 1, fprintf(o.fout,'  No termination message set\n'); end;
if strcmp(r.msg,'opt') == 1, fprintf(o.fout,'  Optimal solution found\n'); end;
if strcmp(r.msg,'inf') == 1, fprintf(o.fout,'  Infeasible stationary point found\n'); end;
if strcmp(r.msg,'itr') == 1, fprintf(o.fout,'  Iteration limit reached\n'); end;
fprintf(o.fout,'\n');

% Print iterate quantities
fprintf(o.fout,'Final values\n');
fprintf(o.fout,'============\n');
fprintf(o.fout,'  Objective function........................ : %+e\n',z.f);
fprintf(o.fout,'  Infeasibility violation................... : %+e\n',z.v);
fprintf(o.fout,'  Optimality error.......................... : %+e\n',z.kkt(2));
fprintf(o.fout,'  Penalty parameter......................... : %+e\n',z.rho);
fprintf(o.fout,'  Interior-point parameter.................. : %+e\n',z.mu);
fprintf(o.fout,'\n');

% Print counters
fprintf(o.fout,'Final counters\n');
fprintf(o.fout,'==============\n');
fprintf(o.fout,'  Iterations................................ : %d\n',c.k);
fprintf(o.fout,'  Function evaluations...................... : %d\n',c.f);
fprintf(o.fout,'  Gradient evaluations...................... : %d\n',c.g);
fprintf(o.fout,'  Hessian evaluations....................... : %d\n',c.H);
fprintf(o.fout,'  Matrix factorizations..................... : %d\n',c.fact);
fprintf(o.fout,'  CPU seconds............................... : %d\n',ceil(toc));
