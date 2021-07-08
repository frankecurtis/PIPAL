% Pipal class
classdef Pipal

  % Class properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    i % Input object
    o % Output object
    c % Counter object
    p % Parameter object
    z % Iterate object
    d % Direction object
    a % Acceptance object
    
  end
  
  % Class methods
  methods
  
    % Constructor
    function P = Pipal(nl,algorithm)

      % Construct classes
      P.i = Input(nl);
      P.o = Output;
      P.c = Counter;
      P.p = Parameter(algorithm);
      P.z = Iterate(P.i,P.c,P.p);
      P.d = Direction;
      P.a = Acceptance;

    end
    
    % Optimization algorithm
    function optimize(P)

      % Print header and line break
      P.o.printHeader(P.i,    P.z);
      P.o.printBreak (    P.c    );
      
      % Iteration loop
      while ~P.z.checkTermination(P.c,P.p)
        P.o.printIterate   (    P.c,    P.z        );
        P.d.evalStep       (P.i,P.c,P.p,P.z,    P.a);
        P.o.printDirection (            P.z,P.d    );
        P.a.lineSearch     (P.i,P.c,P.p,P.z,P.d    );
        P.o.printAcceptance(                    P.a);
        P.z.updateIterate  (P.i,P.c,P.p,    P.d,P.a);
        P.c.incrementIterationCount                 ;
        P.o.printBreak     (    P.c                );
      end
      
      % Print footer and terminate
      P.o.printFooter(P.c,P.p,P.z);
      P.o.terminate               ;
      
    end
        
  end
  
end