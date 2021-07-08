function print_break(o,c)

% function print_break(o,c)
%
% Author       : Frank E. Curtis
% Description  : Prints break in output with header strings.
% Input        : o ~ outputs
%                c ~ counters
% Last revised : 21 June 2010

% Print break in output
if mod(c.k,20) == 0, fprintf(o.fout,'%s\n%s\n%s\n',o.line,o.quan,o.line); end;
