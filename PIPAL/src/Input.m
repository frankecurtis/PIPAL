% Input class
classdef Input < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)

    n0 % number of original formulation variables
    I1 % indices of free variables
    I2 % indices of fixed variables
    I3 % indices of lower bounded variables
    I4 % indices of upper bounded variables
    I5 % indices of lower and upper bounded variables
    I6 % indices of equality constraints
    I7 % indices of lower bounded constraints
    I8 % indices of upper bounded constraints
    I9 % indices of lower and upper bounded constraints
    b2 % right-hand side of fixed variables
    l3 % right-hand side of lower bounded variables
    u4 % right-hand side of upper bounded variables
    l5 % right-hand side of lower half of lower and upper bounded variables
    u5 % right-hand side of upper half of lower and upper bounded variables
    b6 % right-hand side of equality constraints
    l7 % right-hand side of lower bounded constraints
    u8 % right-hand side of upper bounded constraints
    l9 % right-hand side of lower half of lower and upper bounded constraints
    u9 % right-hand side of upper half of lower and upper bounded constraints
    n1 % number of free variables
    n2 % number of fixed variables
    n3 % number of lower bounded variables
    n4 % number of upper bounded variables
    n5 % number of lower and upper bounded variables
    n6 % number of equality constraints
    n7 % number of lower bounded constraints
    n8 % number of upper bounded constraints
    n9 % number of lower and upper bounded constraints
    nV % number of variables
    nI % number of inequality constraints
    nE % number of equality constraints
    nA % size of primal-dual matrix
    x0 % initial point
    
  end
  
  % Class methods
  methods
    
    % Constructor
    function i = Input(nl)
      
      % Initialize AMPL data
      [p.x,p.bl,p.bu,p.l,p.cl,p.cu] = amplfunc(nl);

      % Set number of original formulation variables
      i.n0 = length(p.x);
      
      % Find index sets
      i.I1 = find(p.bl <= -1e18 & p.bu >= 1e18);
      i.I2 = find(p.bl ==  p.bu);
      i.I3 = find(p.bl >  -1e18 & p.bu >= 1e18);
      i.I4 = find(p.bl <= -1e18 & p.bu <  1e18);
      i.I5 = find(p.bl >  -1e18 & p.bu <  1e18 & p.bl ~= p.bu);
      i.I6 = find(p.cl ==  p.cu);
      i.I7 = find(p.cl >  -1e18 & p.cu >= 1e18);
      i.I8 = find(p.cl <= -1e18 & p.cu <  1e18);
      i.I9 = find(p.cl >  -1e18 & p.cu <  1e18 & p.cl ~= p.cu);
      
      % Set right-hand side values
      i.b2 = p.bl(i.I2);
      i.l3 = p.bl(i.I3);
      i.u4 = p.bu(i.I4);
      i.l5 = p.bl(i.I5);
      i.u5 = p.bu(i.I5);
      i.b6 = p.cl(i.I6);
      i.l7 = p.cl(i.I7);
      i.u8 = p.cu(i.I8);
      i.l9 = p.cl(i.I9);
      i.u9 = p.cu(i.I9);
      
      % Set sizes of index sets
      i.n1 = length(i.I1);
      i.n2 = length(i.I2);
      i.n3 = length(i.I3);
      i.n4 = length(i.I4);
      i.n5 = length(i.I5);
      i.n6 = length(i.I6);
      i.n7 = length(i.I7);
      i.n8 = length(i.I8);
      i.n9 = length(i.I9);
      
      % Set number of variables and constraints
      i.nV = i.n1 + i.n3 + i.n4 + i.n5;
      i.nI = i.n3 + i.n4 + 2*i.n5 + i.n7 + i.n8 + 2*i.n9;
      i.nE = i.n6;
      
      % Set size of primal-dual matrix
      i.nA = i.nV + 3*i.nE + 3*i.nI;
      
      % Set initial point
      i.x0 = [p.x(i.I1); p.x(i.I3); p.x(i.I4); p.x(i.I5);];
      
    end
    
  end
  
end
