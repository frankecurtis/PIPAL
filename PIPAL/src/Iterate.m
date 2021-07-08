% Iterate class
classdef Iterate < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    x     % Primal point
    rho   % Penalty parameter value
    rho_  % Penalty parameter last value
    mu    % Interior-point parameter value
    f     % Objective function value (scaled)
    fu    % Objective function value (unscaled)
    g     % Objective gradient value
    r1    % Equality constraint slack value
    r2    % Equality constraint slack value
    cE    % Equality constraint value (scaled)
    JE    % Equality constraint Jacobian value
    JEnnz % Equality constraint Jacobian nonzeros
    lE    % Equality constraint multipliers
    s1    % Inequality constraint slack value
    s2    % Inequality constraint slack value
    cI    % Inequality constraint value (scaled)
    JI    % Inequality constraint Jacobian value
    JInnz % Inequality constraint Jacobian nonzeros
    lI    % Inequality constraint multipliers
    H     % Hessian of Lagrangian
    Hnnz  % Hessian of Lagrangian nonzeros
    v     % Feasibility violation measure value (scaled)
    vu    % Feasibility violation measure value (unscaled)
    v0    % Feasibility violation measure initial value
    phi   % Merit function value
    AL    % Newton matrix L-factor in LDL factorization
    AD    % Newton matrix D-factor in LDL factorization
    AP    % Newton matrix P-factor in LDL factorization
    AS    % Newton matrix S-factor in LDL factorization
    shift % Hessian shift value
    b     % Newton right-hand side
    kkt   % KKT errors
    kkt_  % KKT errors last value
    
  end
  
  % Class properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    fs      % Objective scaling factor
    cEs     % Equality constraint scaling factors
    cEu     % Equality constraint value (unscaled)
    cIs     % Inequality constraint scaling factors
    cIu     % Inequality constraint value (unscaled)
    A       % Newton matrix
    shift22 % Newton matrix (2,2)-block shift value
    v_      % Feasibility violation measure last value
    cut_    % Boolean value for last backtracking line search

  end
  
  % Class methods
  methods

    % Constructor
    function z = Iterate(i,c,p)
      
      % Initialize quantities
      z.x                 = i.x0;
      z.rho               = p.rho_init;
      z.mu                = p.mu_init;
      z.lE                = zeros(i.nE,1);
      z.lI                = (1/2)*ones(i.nI,1);
      z.evalScalings        (i,c,p);
      z.evalFunctions       (i,c  );
      z.evalGradients       (i,c  );
      z.evalDependent       (i,  p);
      z.v0                = 1;
      z.evalInfeasibility   (i    );
      z.v0                = z.v;
      z.evalInfeasibility   (i    );
      z.v_                = z.v;
      z.shift             = 0;
      z.kkt_              = inf*ones(p.opt_err_mem,1);
      z.cut_              = 0;
      z.evalHessian         (i,c  );
      z.Hnnz              = nnz(z.H);
      z.JEnnz             = nnz(z.JE);
      z.JInnz             = nnz(z.JI);
      z.initNewtonMatrix    (i    );
      z.evalNewtonMatrix    (i,c,p);
      
    end
    
    % Termination checker
    %   References : c.*, p.*, kkt(1:2), v, v0
    %   Calls      : -
    %   Modifies   : -
    function b = checkTermination(z,c,p)

      % Initialize boolean
      b = 0;
      
      % Update termination based on optimality error of nonlinear optimization problem
      if z.kkt(2) <= p.opt_err_tol & z.v <= p.opt_err_tol, b = 1; return; end;
      
      % Update termination based on optimality error of feasibility problem
      if z.kkt(1) <= p.opt_err_tol & z.v >  p.opt_err_tol, b = 2; return; end;
      
      % Update termination based on iteration count
      if c.k >= p.iter_max, b = 3; return; end;
      
    end
    
    % Dependent quantity evaluator
    %   References : i.*, p.*, rho, mu, f, g, cE, JE, r1, r2, lE, cI, JI, s1, s2, lI
    %   Calls      : evalSlacks, evalMerit, evalKKTErrors, evalKKTError
    %   Modifies   : r1, r2, s1, s2, phi, kkt
    function     evalDependent(z,i,p)
      
      % Evaluate quantities dependent on penalty and interior-point parameters
      z.evalSlacks   (i,p);
      z.evalMerit    (i  );
      z.evalKKTErrors(i  );

    end
    
    % Function evaluator
    %   References : i.*, x, fs, cEs, cIs
    %   Calls      : evalXOriginal, c.incrementFunctionCount, amplfunc
    %   Modifies   : f, cE, cI, fu, cEu, cIu
    function     evalFunctions(z,i,c)
      
      % Evaluate x in original space
      x_orig = z.evalXOriginal(i);
      
      % Evaluate AMPL functions
      [z.f,c_orig] = amplfunc(x_orig,0);

      % Increment function evaluation counter
      c.incrementFunctionCount;
      
      % Set equality constraint values
      if i.nE > 0, z.cE = c_orig(i.I6) - i.b6; end;
      
      % Initialize inequality constraint values
      if i.nI > 0, z.cI = zeros(i.nI,1); end;
      
      % Set inequality constraint values
      if i.n3 > 0, z.cI(1                                   :i.n3                                   ) =  i.l3 - z.x(1+i.n1          :i.n1+i.n3          ); end;
      if i.n4 > 0, z.cI(1+i.n3                              :i.n3+i.n4                              ) = -i.u4 + z.x(1+i.n1+i.n3     :i.n1+i.n3+i.n4     ); end;
      if i.n5 > 0, z.cI(1+i.n3+i.n4                         :i.n3+i.n4+i.n5                         ) =  i.l5 - z.x(1+i.n1+i.n3+i.n4:i.n1+i.n3+i.n4+i.n5);
                   z.cI(1+i.n3+i.n4+i.n5                    :i.n3+i.n4+i.n5+i.n5                    ) = -i.u5 + z.x(1+i.n1+i.n3+i.n4:i.n1+i.n3+i.n4+i.n5); end;
      if i.n7 > 0, z.cI(1+i.n3+i.n4+i.n5+i.n5               :i.n3+i.n4+i.n5+i.n5+i.n7               ) =  i.l7 - c_orig(i.I7);                              end;
      if i.n8 > 0, z.cI(1+i.n3+i.n4+i.n5+i.n5+i.n7          :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8          ) = -i.u8 + c_orig(i.I8);                              end;
      if i.n9 > 0, z.cI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8     :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9     ) =  i.l9 - c_orig(i.I9);
                   z.cI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9:i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9) = -i.u9 + c_orig(i.I9);                              end;
      
      % Store unscaled quantities
      z.fu = z.f;
      if i.nE > 0, z.cEu = z.cE; end;
      if i.nI > 0, z.cIu = z.cI; end;
      
      % Scale quantities
      z.f = z.fs*z.f;
      if i.nE > 0, z.cE = z.cEs.*z.cE; end;
      if i.nI > 0, z.cI = z.cIs.*z.cI; end;
      
    end
    
    % Gradient evaluator
    %   References : i.*, x, fs, cEs, cIs
    %   Calls      : evalXOriginal, c.incrementGradientCount, amplfunc
    %   Modifies   : g, JE, JI
    function     evalGradients(z,i,c)
      
      % Evaluate x in original space
      x_orig = z.evalXOriginal(i);
      
      % Evaluate AMPL gradients
      [g_orig,J_orig] = amplfunc(x_orig,1);

      % Increment gradient evaluation counter
      c.incrementGradientCount;
      
      % Set objective gradient
      z.g = [g_orig(i.I1); g_orig(i.I3); g_orig(i.I4); g_orig(i.I5)];
      
      % Set equality constraint Jacobian
      if i.nE > 0, z.JE = [J_orig(i.I6,i.I1) J_orig(i.I6,i.I3) J_orig(i.I6,i.I4) J_orig(i.I6,i.I5)]; end;
      
      % Initialize inequality constraint Jacobian
      if i.nI > 0, z.JI = sparse(i.nI,i.nV); end;
      
      % Set inequality constraint Jacobian
      if i.n3 > 0, z.JI(1                                   :i.n3                                   ,1+i.n1          :i.n1+i.n3          ) = -speye(i.n3); end;
      if i.n4 > 0, z.JI(1+i.n3                              :i.n3+i.n4                              ,1+i.n1+i.n3     :i.n1+i.n3+i.n4     ) =  speye(i.n4); end;
      if i.n5 > 0, z.JI(1+i.n3+i.n4                         :i.n3+i.n4+i.n5                         ,1+i.n1+i.n3+i.n4:i.n1+i.n3+i.n4+i.n5) = -speye(i.n5);
                   z.JI(1+i.n3+i.n4+i.n5                    :i.n3+i.n4+i.n5+i.n5                    ,1+i.n1+i.n3+i.n4:i.n1+i.n3+i.n4+i.n5) =  speye(i.n5); end;
      if i.n7 > 0, z.JI(1+i.n3+i.n4+i.n5+i.n5               :i.n3+i.n4+i.n5+i.n5+i.n7               ,1               :i.n1+i.n3+i.n4+i.n5) = -[J_orig(i.I7,i.I1) J_orig(i.I7,i.I3) J_orig(i.I7,i.I4) J_orig(i.I7,i.I5)]; end;
      if i.n8 > 0, z.JI(1+i.n3+i.n4+i.n5+i.n5+i.n7          :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8          ,1               :i.n1+i.n3+i.n4+i.n5) =  [J_orig(i.I8,i.I1) J_orig(i.I8,i.I3) J_orig(i.I8,i.I4) J_orig(i.I8,i.I5)]; end;
      if i.n9 > 0, z.JI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8     :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9     ,1               :i.n1+i.n3+i.n4+i.n5) = -[J_orig(i.I9,i.I1) J_orig(i.I9,i.I3) J_orig(i.I9,i.I4) J_orig(i.I9,i.I5)];
                   z.JI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9:i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9,1               :i.n1+i.n3+i.n4+i.n5) =  [J_orig(i.I9,i.I1) J_orig(i.I9,i.I3) J_orig(i.I9,i.I4) J_orig(i.I9,i.I5)]; end;
      
      % Scale objective gradient
      z.g = z.fs*z.g;
      
      % Scale constraint Jacobians
      if i.nE > 0, z.JE = spdiags(z.cEs,0:0,i.nE,i.nE)*z.JE; end;
      if i.nI > 0, z.JI = spdiags(z.cIs,0:0,i.nI,i.nI)*z.JI; end;

    end
    
    % Hessian evaluator
    %   References : i.*, x, fs, cEs, cIs, lE, lI, rho
    %   Calls      : evalXOriginal, c.incrementHessianCount, amplfunc
    %   Modifies   : H
    function     evalHessian(z,i,c)

      % Evaluate x in original space
      x_orig = z.evalXOriginal(i);
      
      % Initialize multipliers in original space
      l_orig = zeros(i.nE+i.n7+i.n8+i.n9,1);
      
      % Scale equality constraint multipliers
      if i.nE > 0, lE = z.lE.*(z.cEs/(z.rho*z.fs)); end;
      
      % Set equality constraint multipliers in original space
      if i.nE > 0, l_orig(i.I6) = lE; end;
      
      % Scale inequality constraint multipliers
      if i.n7+i.n8+i.n9 > 0, lI = z.lI.*(z.cIs/(z.rho*z.fs)); end;
      
      % Set inequality constraint multipliers in original space
      if i.n7 > 0, l_orig(i.I7) = -lI(1+i.n3+i.n4+i.n5+i.n5               :i.n3+i.n4+i.n5+i.n5+i.n7               ); end;
      if i.n8 > 0, l_orig(i.I8) = +lI(1+i.n3+i.n4+i.n5+i.n5+i.n7          :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8          ); end;
      if i.n9 > 0, l_orig(i.I9) = -lI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8     :i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9     )...
                                  +lI(1+i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9:i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9); end;
      
      % Evaluate H_orig
      if (i.nE+i.n7+i.n8+i.n9 == 0), H_orig = amplfunc([]);
      else                           H_orig = amplfunc(l_orig); end;
      
      % Increment Hessian evaluation counter
      c.incrementHessianCount;
      
      % Set Hessian of the Lagrangian
      z.H = [H_orig(i.I1,i.I1) H_orig(i.I1,i.I3) H_orig(i.I1,i.I4) H_orig(i.I1,i.I5);
             H_orig(i.I3,i.I1) H_orig(i.I3,i.I3) H_orig(i.I3,i.I4) H_orig(i.I3,i.I5);
             H_orig(i.I4,i.I1) H_orig(i.I4,i.I3) H_orig(i.I4,i.I4) H_orig(i.I4,i.I5);
             H_orig(i.I5,i.I1) H_orig(i.I5,i.I3) H_orig(i.I5,i.I4) H_orig(i.I5,i.I5);];
      
      % Rescale H and Ht
      z.H = z.rho*z.fs*z.H;
      
    end
    
    % Infeasibility evaluator
    %   References : i.*, cE, cI, v0
    %   Calls      : evalViolation
    %   Modifies   : v, vu
    function     evalInfeasibility(z,i)
      
      % Evaluate scaled and unscaled feasibility violations
      z.v  = z.evalViolation(i,z.cE ,z.cI )/max(1,z.v0);
      z.vu = z.evalViolation(i,z.cEu,z.cIu);
      
    end
    
    % KKT error evaluator
    %   References : i.*, g, JE, lE, JI, lI, r1, r2, s1, s2
    %   Calls      : -
    %   Modifies   : -
    function v = evalKKTError(z,i,rho,mu)
      
      % Initialize optimality vector
      kkt = zeros(i.nV+2*i.nE+2*i.nI,1);
      
      % Set gradient of penalty objective
      kkt(1:i.nV) = rho*z.g;
      
      % Set gradient of Lagrangian for constraints
      if i.nE > 0, kkt(1:i.nV) = kkt(1:i.nV) + (z.lE'*z.JE)'; end;
      if i.nI > 0, kkt(1:i.nV) = kkt(1:i.nV) + (z.lI'*z.JI)'; end;
      
      % Set complementarity for constraint slacks
      if i.nE > 0, kkt(1+i.nV       :i.nV+2*i.nE       ) = [z.r1.*(1 + z.lE) - mu; z.r2.*(1 - z.lE) - mu]; end;
      if i.nI > 0, kkt(1+i.nV+2*i.nE:i.nV+2*i.nE+2*i.nI) = [z.s1.*(0 + z.lI) - mu; z.s2.*(1 - z.lI) - mu]; end;
      
      % Scale complementarity
      if rho > 0, kkt = (1/max(1,norm(rho*z.g,inf)))*kkt; end;
      
      % Evaluate optimality error
      v = norm(kkt,inf);

    end
    
    % KKT errors evaluator
    %   References : i.*, rho, mu, g, JE, lE, JI, lI, r1, r2, s1, s2
    %   Calls      : evalKKTError
    %   Modifies   : kkt(1:3)
    function     evalKKTErrors(z,i)
      
      % Loop to compute optimality errors
      z.kkt(1) = z.evalKKTError(i,0    ,0   );
      z.kkt(2) = z.evalKKTError(i,z.rho,0   );
      z.kkt(3) = z.evalKKTError(i,z.rho,z.mu);
      
    end

    % Matrices evaluator
    %   References : i.*, p.*, x, fs, cEs, cIs, rho, H, JE, r1, r2, lE, JI, s1, s2, lI, shift, shift22, cut_
    %   Calls      : evalHessian, evalXOriginal, c.incrementHessianCount, amplfunc, evalNewtonMatrix, c.incrementFactorizationCount
    %   Modifies   : H, A, AL, AD, AP, AS, shift, shift22
    function     evalMatrices(z,i,c,p)

      % Evaluate Hessian and Newton matrices
      z.evalHessian     (i,c  );
      z.evalNewtonMatrix(i,c,p);

    end
    
    % Merit evaluator
    %   References : i.*, rho, f, mu, r1, r2, s1, s2
    %   Calls      : -
    %   Modifies   : phi
    function     evalMerit(z,i)
      
      % Initialize merit for objective
      z.phi = z.rho*z.f;

      % Update merit for slacks
      if i.nE > 0, z.phi = z.phi - z.mu*sum(log([z.r1;z.r2])) + sum([z.r1;z.r2]); end;
      if i.nI > 0, z.phi = z.phi - z.mu*sum(log([z.s1;z.s2])) + sum(      z.s2 ); end;
      
    end
    
    % Newton matrix evaluator
    %   References : i.*, p.*, H, JE, r1, r2, lE, JI, s1, s2, lI, shift, shift22, cut_
    %   Calls      : c.incrementFactorizationCount
    %   Modifies   : H, A, AL, AD, AP, AS, shift, shift22
    function     evalNewtonMatrix(z,i,c,p)
      
      % Check for equality constraints
      if i.nE > 0
        
        % Set diagonal terms
        for j = 1:i.nE, z.A(i.nV+     j,i.nV+     j) = (1+z.lE(j))/z.r1(j);
                        z.A(i.nV+i.nE+j,i.nV+i.nE+j) = (1-z.lE(j))/z.r2(j); end;
          
        % Set constraint Jacobian
        z.A(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI,1:i.nV) = z.JE;
        
      end
      
      % Check for inequality constraints
      if i.nI > 0
          
        % Set diagonal terms
        for j = 1:i.nI, z.A(i.nV+2*i.nE+     j,i.nV+2*i.nE+     j) = (0+z.lI(j))/z.s1(j);
                        z.A(i.nV+2*i.nE+i.nI+j,i.nV+2*i.nE+i.nI+j) = (1-z.lI(j))/z.s2(j); end;
          
        % Set constraint Jacobian
        z.A(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI,1:i.nV) = z.JI;

      end
      
      % Set minimum potential shift
      min_shift = max(p.shift_min,p.shift_factor*z.shift);
      
      % Initialize Hessian modification
      if z.cut_ == 1, z.shift = min(p.shift_max,min_shift/p.shift_factor); else z.shift = 0; end;
      
      % Initialize inertia correction loop
      done = 0; z.shift22 = 0;
      
      % Loop until inertia is correct
      while ~done & z.shift < p.shift_max
      
        % Set Hessian of Lagrangian
        z.A(1:i.nV,1:i.nV) = z.H+z.shift*speye(i.nV);
        
        % Set diagonal terms
        for j = 1:i.nE, z.A(i.nV+2*i.nE+2*i.nI+j,i.nV+2*i.nE+2*i.nI+j) = -z.shift22; end;
        
        % Set diagonal terms
        for j = 1:i.nI, z.A(i.nV+3*i.nE+2*i.nI+j,i.nV+3*i.nE+2*i.nI+j) = -z.shift22; end;
        
        % Factor primal-dual matrix
        [z.AL,z.AD,z.AP,z.AS,neig] = ldl(tril(z.A),p.pivot_thresh,'vector');
        
        % Increment factorization counter
        c.incrementFactorizationCount;
        
        % Set number of nonnegative eigenvalues
        peig = i.nA - neig;
        
        % Check inertia
        if     peig < i.nV+2*i.nE+2*i.nI        , z.shift   = max(min_shift,z.shift/p.shift_factor);
        elseif neig < i.nE+i.nI & z.shift22 == 0, z.shift22 = p.shift_min;
        else                                      done      = 1; end;
        
      end
      
      % Update Hessian
      z.H = z.H+z.shift*speye(i.nV);
      
    end
    
    % Newton right-hand side evaluator
    %   References : i.*, rho, mu, g, cE, JE, r1, r2, lE, cI, JI, s1, s2, lI
    %   Calls      : -
    %   Modifies   : b
    function     evalNewtonRhs(z,i)
      
      % Initialize right-hand side vector
      z.b = zeros(i.nA,1);
      
      % Set gradient of objective
      z.b(1:i.nV) = z.rho*z.g;
      
      % Set gradient of Lagrangian for constraints
      if i.nE > 0, z.b(1:i.nV) = z.b(1:i.nV) + (z.lE'*z.JE)'; end;
      if i.nI > 0, z.b(1:i.nV) = z.b(1:i.nV) + (z.lI'*z.JI)'; end;
      
      % Set complementarity for constraint slacks
      if i.nE > 0, z.b(1+i.nV       :i.nV+2*i.nE       ) = [1 + z.lE - z.mu./z.r1; 1 - z.lE - z.mu./z.r2]; end;
      if i.nI > 0, z.b(1+i.nV+2*i.nE:i.nV+2*i.nE+2*i.nI) = [0 + z.lI - z.mu./z.s1; 1 - z.lI - z.mu./z.s2]; end;
      
      % Set penalty-interior-point constraint values
      if i.nE > 0, z.b(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI) = z.cE + z.r1 - z.r2; end;
      if i.nI > 0, z.b(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI) = z.cI + z.s1 - z.s2; end;
      
    end
    
    % Scalings evaluator
    %   References : i.*, p.*, x, g, JE, JI, fs, cEs, cIs
    %   Calls      : evalGradients, evalXOriginal, c.incrementGradientCount, amplfunc
    %   Modifies   : g, JE, JI, fs, cEs, cIs
    function     evalScalings(z,i,c,p)
      
      % Initialize scalings
      z.fs  = 0;
      z.cEs = zeros(i.nE,1);
      z.cIs = zeros(i.nI,1);
      
      % Evaluate gradients
      z.evalGradients(i,c);

      % Evaluate norm of objective gradient
      g_norm_inf = norm(z.g,inf);
      
      % Scale down objective if norm of gradient is too large
      z.fs = p.grad_max/max(g_norm_inf,p.grad_max);
      
      % Loop through equality constraints
      for j = 1:i.nE
      
        % Evaluate norm of gradient of jth equality constraint
        JE_j_norm_inf = norm(z.JE(j,:),inf);
        
        % Scale down equality constraint j if norm of gradient is too large
        z.cEs(j) = p.grad_max/max(JE_j_norm_inf,p.grad_max);
      
      end
      
      % Loop through inequality constraints
      for j = 1:i.nI
      
        % Evaluate norm of gradient of jth inequality constraint
        JI_j_norm_inf = norm(z.JI(j,:),inf);
        
        % Scale down inequality constraint j if norm of gradient is too large
        z.cIs(j) = p.grad_max/max(JI_j_norm_inf,p.grad_max);
      
      end
      
    end
    
    % Slacks evaluator
    %   References : i.*, p.*, mu, cE, cI
    %   Calls      : -
    %   Modifies   : r1, r2, s1, s2
    function     evalSlacks(z,i,p)
      
      % Check for equality constraints
      if i.nE > 0

        % Set slacks
        z.r1 = (1/2)*(z.mu - z.cE + sqrt(z.cE.^2 + z.mu^2));
        z.r2 = (1/2)*(z.mu + z.cE + sqrt(z.cE.^2 + z.mu^2));
        
        % Adjust for numerical error
        z.r1 = max(z.r1,p.slack_min);
        z.r2 = max(z.r2,p.slack_min);
        
      end
      
      % Check for inequality constraints
      if i.nI > 0
        
        % Set slacks
        z.s1 = (1/2)*(2*z.mu - z.cI + sqrt(z.cI.^2 + 4*z.mu^2));
        z.s2 = (1/2)*(2*z.mu + z.cI + sqrt(z.cI.^2 + 4*z.mu^2));

        % Adjust for numerical error
        z.s1 = max(z.s1,p.slack_min);
        z.s2 = max(z.s2,p.slack_min);
        
      end
      
    end
    
    % Evaluator of x in original space
    %   References : i.*, x
    %   Calls      : -
    %   Modifies   : -
    function x = evalXOriginal(z,i)
      
      % Initialize x in original space
      x = zeros(i.n0,1);
      
      % Evaluate x in original space
      x(i.I1) = z.x(1               :i.n1               );
      x(i.I2) = i.b2;
      x(i.I3) = z.x(1+i.n1          :i.n1+i.n3          );
      x(i.I4) = z.x(1+i.n1+i.n3     :i.n1+i.n3+i.n4     );
      x(i.I5) = z.x(1+i.n1+i.n3+i.n4:i.n1+i.n3+i.n4+i.n5);
      
    end

    % Initializes Newton matrix
    %   References : i.*, Hnnz, JEnnz, JInnz
    %   Calls      : -
    %   Modifies   : A
    function     initNewtonMatrix(z,i)

      % Allocate memory
      z.A = spalloc(i.nA,i.nA,z.Hnnz+5*i.nE+5*i.nI+z.JEnnz+z.JInnz);

      % Initialize interior-point Hessians
      z.A(1+i.nV       :i.nV+2*i.nE       ,1+i.nV       :i.nV+2*i.nE       ) = speye(2*i.nE);
      z.A(1+i.nV+2*i.nE:i.nV+2*i.nE+2*i.nI,1+i.nV+2*i.nE:i.nV+2*i.nE+2*i.nI) = speye(2*i.nI);

      % Check for equality constraints
      if i.nE > 0
                  
        % Initialize constraint Jacobian
        z.A(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI,1+i.nV              :i.nV+  i.nE       ) =  speye(i.nE);
        z.A(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI,1+i.nV+  i.nE       :i.nV+2*i.nE       ) = -speye(i.nE);
        z.A(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI,1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI) = -speye(i.nE);
        
      end
        
      % Check for inequality constraints
      if i.nI > 0
          
        % Initialize constraint Jacobian
        z.A(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI,1+i.nV+2*i.nE       :i.nV+2*i.nE+  i.nI) =  speye(i.nI);
        z.A(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI,1+i.nV+2*i.nE+  i.nI:i.nV+2*i.nE+2*i.nI) = -speye(i.nI);
        z.A(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI,1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI) = -speye(i.nI);
        
      end

    end

    % Set interior-point parameter
    %   References : -
    %   Calls      : -
    %   Modifies   : mu
    function     setMu(z,mu)
      z.mu = mu;
    end

    % Set primal variables
    %   References : -
    %   Calls      : -
    %   Modifies   : x, r1, r2, lE, s1, s2, lI
    function     setPrimals(z,i,x,r1,r2,s1,s2,lE,lI)

      % Set primal variables
      z.x = x;
      if i.nE > 0, z.r1 = r1; z.r2 = r2; z.lE = lE; end;
      if i.nI > 0, z.s1 = s1; z.s2 = s2; z.lI = lI; end;

    end
    
    % Set penalty parameter
    %   References : -
    %   Calls      : -
    %   Modifies   : rho
    function     setRho(z,rho)
      z.rho = rho;
    end
    
    % Set last penalty parameter
    %   References : -
    %   Calls      : -
    %   Modifies   : rho_
    function     setRhoLast(z,rho)
      z.rho_ = rho;
    end
    
    % Iterate updater
    %   References : i.*, p.*, a.p, a.p0, a.d, d.x, d.r1, d.r2, d.lE, d.s1, d.s2, d.lI, rho, mu, x, f, fs, v, v0, g, cE, cEs, JE, r1, r2, lE, cI, cIs, JI, s1, s2, lI, kkt
    %   Calls      : updatePoint, evalInfeasibility, evalViolation, evalGradients, evalXOriginal, c.incrementGradientCount, amplfunc, evalSlacks, evalMerit, evalKKTErrors, evalKKTError
    %   Modifies   : v_, cut_, kkt_, x, v, vu, g, JE, r1, r2, lE, JI, s1, s2, lI, phi, kkt
    function     updateIterate(z,i,c,p,d,a)
      
      % Update last quantities
      z.v_   = z.v;
      z.cut_ = (a.p < a.p0);
      
      % Update iterate quantities
      z.updatePoint      (i,    d,a);
      z.evalInfeasibility(i        );
      z.evalGradients    (i,c      );
      z.evalDependent    (i,  p    );
      
      % Update last KKT errors
      z.kkt_ = [z.kkt(2); z.kkt_(1:p.opt_err_mem-1)];

    end
    
    % Parameter updater
    %   References : i.*, p.*, rho, mu, f, v, v_, g, cE, JE, r1, r2, lE, cI, JI, s1, s2, lI, kkt
    %   Calls      : p.setMuMaxExpZero, evalDependent, evalSlacks, evalMerit, evalKKTErrors, evalKKTError
    %   Modifies   : p.mu_max_exp, rho, mu, r1, r2, s1, s2, phi, kkt
    function     updateParameters(z,i,p)

      % Check for interior-point parameter update based on optimality error
      while z.mu > p.mu_min & z.kkt(3) <= max([z.mu;p.opt_err_tol-z.mu])
      
        % Restrict interior-point parameter increase
        p.setMuMaxExpZero;
        
        % Update interior-point parameter
        if z.mu > p.mu_min
        
          % Decrease interior-point
          z.mu = max(p.mu_min,p.mu_factor*z.mu);
          
          % Evaluate penalty and interior-point parameter dependent quantities
          z.evalDependent(i,p);

        end
        
      end
      
      % Check for penalty parameter update based on optimality error
      if (z.kkt(2) <= p.opt_err_tol & z.v > p.opt_err_tol) | z.v > max([1;z.v_;p.infeas_max])
        
        % Update penalty parameter
        if z.rho > p.rho_min
        
          % Decrease penalty parameter
          z.rho = max(p.rho_min,p.rho_factor*z.rho);
          
          % Evaluate penalty and interior-point parameter dependent quantities
          z.evalDependent(i,p);
        
        end
        
      end

    end
    
    % Primal point updater
    %   References : i.*, x, r1, r2, lE, s1, s2, lI, a.p, a.d, d.x, d.r1, d.r2, d.lE, d.s1, d.s2, d.lI
    %   Calls      : -
    %   Modifies   : x, r1, r2, lE, s1, s2, lI
    function     updatePoint(z,i,d,a)
      
      % Update primal and dual variables
                   z.x  = z.x  + a.p*d.x ;
      if i.nE > 0, z.r1 = z.r1 + a.p*d.r1;
                   z.r2 = z.r2 + a.p*d.r2; end;
      if i.nI > 0, z.s1 = z.s1 + a.p*d.s1;
                   z.s2 = z.s2 + a.p*d.s2; end;
      if i.nE > 0, z.lE = z.lE + a.d*d.lE; end;
      if i.nI > 0, z.lI = z.lI + a.d*d.lI; end;
      
    end
    
  end
  
  % Class methods (static)
  methods (Static)

    % Feasibility violation evaluator
    %   References : i.*
    %   Calls      : -
    %   Modifies   : -
    function v = evalViolation(i,cE,cI)
      
      % Initialize violation vector
      vec = [];
      
      % Update vector for constraint values
      if i.nE > 0, vec = cE;               end;
      if i.nI > 0, vec = [vec; max(cI,0)]; end;
      
      % Evaluate vector norm
      v = norm(vec,1);

    end
    
  end
  
end
