% Acceptance class
classdef Acceptance < handle

  % Class properties (private set access)
  properties (SetAccess = private)

    p0 = 0; % Fraction-to-the-boundary steplength
    p  = 0; % Primal steplength
    d  = 0; % Dual steplength
    
  end
  
  % Class methods
  methods
    
    % Backtracking line search
    function backtracking(a,i,c,p,z,d)
      
      % Store current values
      phi = z.phi; x = z.x;
      r1 = z.r1; r2 = z.r2; lE = z.lE;
      s1 = z.s1; s2 = z.s2; lI = z.lI;
      
      % Backtracking loop
      while a.p >= eps
        
        % Set trial point
        z.updatePoint  (i,    d,a);
        z.evalFunctions(i,c      );
        z.evalSlacks   (i,  p    );
        z.evalMerit    (i        );

        % Check for nonlinear fraction-to-boundary violation
        ftb = 0;
        if i.nE > 0, ftb = ftb + sum(z.r1<p.ls_frac*r1) + sum(z.r2<p.ls_frac*r2); end;
        if i.nI > 0, ftb = ftb + sum(z.s1<p.ls_frac*s1) + sum(z.s2<p.ls_frac*s2); end;
        
        % Reset variables
        z.setPrimals(i,x,r1,r2,s1,s2,lE,lI);
        
        % Check Armijo condition
        if ftb == 0 & z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0), return; end;
        
        % Reduce steplength
        a.p = p.ls_factor*a.p;
        
      end
      
    end
    
    % Fraction-to-boundary line search
    function fractionToBoundary(a,i,p,z,d)
      
      % Initialize primal fraction-to-boundary
      a.p0 = 1;
      
      % Update primal fraction-to-boundary for constraint slacks
      if i.nE > 0, a.p0 = min([a.p0; (p.ls_frac-1)*z.r1(d.r1<0)./d.r1(d.r1<0); (p.ls_frac-1)*z.r2(d.r2<0)./d.r2(d.r2<0)]); end;
      if i.nI > 0, a.p0 = min([a.p0; (p.ls_frac-1)*z.s1(d.s1<0)./d.s1(d.s1<0); (p.ls_frac-1)*z.s2(d.s2<0)./d.s2(d.s2<0)]); end;
      
      % Initialize primal steplength
      a.p = a.p0;
      
      % Initialize dual fraction-to-boundary
      a.d = 1;
      
      % Update dual fraction-to-boundary for constraint multipliers
      if i.nE > 0, a.d = min([a.d; (p.ls_frac-1)*(1+z.lE(d.lE<0))./d.lE(d.lE<0); (1-p.ls_frac)*(1-z.lE(d.lE>0))./d.lE(d.lE>0)]); end;
      if i.nI > 0, a.d = min([a.d; (p.ls_frac-1)*(0+z.lI(d.lI<0))./d.lI(d.lI<0); (1-p.ls_frac)*(1-z.lI(d.lI>0))./d.lI(d.lI>0)]); end;
      
    end
    
    % Full step search for trial penalty parameters
    function b = fullStepCheck(a,i,c,p,z,d)
      
      % Initialize boolean
      b = 0;
      
      % Set current and last penalty parameters
      rho      = z.rho;
      rho_temp = z.rho_;
      
      % Loop through last penalty parameters
      while rho < rho_temp
        
        % Set penalty parameter
        z.setRho(rho_temp);
        
        % Evaluate merit
        z.evalMerit(i);
        
        % Store current values
        phi = z.phi; x = z.x;
        r1 = z.r1; r2 = z.r2; lE = z.lE;
        s1 = z.s1; s2 = z.s2; lI = z.lI;
        
        % Set trial point
        z.updatePoint  (i,    d,a);
        z.evalFunctions(i,c      );
        z.evalSlacks   (i,  p    );
        z.evalMerit    (i        );
        
        % Check for nonlinear fraction-to-boundary violation
        ftb = 0;
        if i.nE > 0, ftb = ftb + sum(z.r1<p.ls_frac*r1) + sum(z.r2<p.ls_frac*r2); end;
        if i.nI > 0, ftb = ftb + sum(z.s1<p.ls_frac*s1) + sum(z.s2<p.ls_frac*s2); end;
        
        % Reset variables
        z.setPrimals(i,x,r1,r2,s1,s2,lE,lI);
        
        % Check Armijo condition
        if ftb == 0 & z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0), b = 1; return; end;
        
        % Decrease rho
        rho_temp = p.rho_factor*rho_temp;
        
      end
      
      % Set rho
      z.setRho(rho);
      
      % Evaluate merit
      z.evalMerit(i);
      
    end
    
    % Line search
    function lineSearch(a,i,c,p,z,d)
      
      % Check fraction-to-boundary rule
      a.fractionToBoundary(i,p,z,d);
      
      % Check for full step for trial penalty parameters
      b = a.fullStepCheck(i,c,p,z,d);
      
      % Run backtracking line search
      if b == 0, a.backtracking(i,c,p,z,d); end;
      
    end
    
  end
  
end