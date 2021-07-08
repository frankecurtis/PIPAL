function o = init_output(i)

% function o = init_output(i)
%
% Author       : Frank E. Curtis
% Description  : Initializes output by opening output file
%                and storing output header strings.
% Input        : i ~ inputs
% Output       : o ~ outputs
% Last revised : 21 June 2010

% Check for output file
if ~isfield(i,'output_file'), i.output_file = 'pipal.out'; end;

% Open output file for writing
o.fout = fopen(i.output_file,'w');

% Error check for writable output file
assert(o.fout~=-1, sprintf('PIPAL: Error opening output file, %s.',i.output_file));

% Check for warnings option
if ~isfield(i,'warnings'), o.warnings = 0; else o.warnings = i.warnings; end;

% Check warnings options
if o.warnings == 0, warning off all; end;

% Set verbosity level
if ~isfield(i,'verbosity'), o.verbosity = 0; else o.verbosity = i.verbosity; end;

% Check verbosity level
if o.verbosity <= 0

  % Store output strings
  o.line = '======+=========================+====================================+============+===========';
  o.quan = 'Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error | ||P.Step|| | Pri. Step.';
  o.none =                                                                        '---------- | ----------';
  
else

  % Store output strings
  o.line = '======+=========================+====================================+=========================+===========================================================================+=======================';
  o.quan = 'Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error |    Merit     P.I.P. Err.|    Shift    ||P.Step||  ||D.Step||   Lin. Red.    Quad. Red.    Quality   | Pri. Step.  Dual Step.';
  o.none =                                                                        '-----------  ---------- | ----------  ----------  ----------  -----------  -----------  ----------- | ----------  ----------';

end
