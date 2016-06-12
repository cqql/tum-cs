% ================================================================
% Scheme and matrices follow Michael Dumbser's ADER DG Code.
%
% Evolves PDE in quasi-linear form
% u_t + A * q_x = 0
%
% ================================================================
%
close all;
clear all;

% f(x>=0) = 1; f(x<0) = 0;
heaviside=@(x)(ceil(sign(x)/2));

global N;

%%%%%%%%%%%%%%%%%%%%%%
% CHANGE IF YOU WANT %
%%%%%%%%%%%%%%%%%%%%%%
N=6;                                        % order of spatial discretization
cells=20;                                   % number of elements
tEnd  = 1.0;                                % end of simulation time
A = [1];                                    % speed at which wave travels (must be > 0)

%%%%%%%%%%%%%%%%%%%%%%
% DO NOT CHANGE      %
%%%%%%%%%%%%%%%%%%%%%%
CFLmax = 0.5 / (2*N + 1);
CFL = 0.9;                                  % Relative Courant number (w.r.t. chosen DG scheme)
IMAX  = cells;                              % number of elements
ghosts = 1;                                 % periodic boundaries (we don't use reconstruction)
nDOFs = N+1;                                % number of spatial degrees of freedom
nDOF  = (N+1)^2;                            % number of space-time degrees of freedom
lmax = max(A(1,1));                         % maximum eigenvalue - trivial for our 1x1 matrix
x0    = 0;                                  % start of the computational domain
x1    = 1;                                  % end of the computational domain
dx    = (x1-x0)/IMAX;                       % mesh spacing (here equispaced)
x     = linspace(x0+dx/2,x1-dx/2,IMAX);


%%%%%%%%%%%%%%%%%%%%%%
% nodes and weights  %
%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Legendre quadrature nodes and weights
% Gaussian quadrature is accurate up to order 2n+1
[xGPN,wGPN] = gaussLegendre(N+1,0,1);
nGP  = length(xGPN);                        % number of quadrature points

gaussPointsGlobal = zeros(cells*nGP,1);
cnt = 1;
size(gaussPointsGlobal(cnt:cnt+nGP-1))
for i=1:IMAX
  gaussPointsGlobal(cnt:cnt+nGP-1) = (i-1)*dx+xGPN*dx;
  cnt = cnt+nGP;
end

%%%%%%%%%%%%%%%%%%%%%%
% data structures    %
%%%%%%%%%%%%%%%%%%%%%%
exact_solution = zeros(IMAX,nGP);
uhat = zeros(nDOFs,IMAX+ghosts*2);          % the DG polynomial

%%%%%%%%%%%%%%%%%%%%%%
% System matrices    %
%%%%%%%%%%%%%%%%%%%%%%
M     = zeros(nDOF,nDOF);                   % mass matrix
Kxi   = zeros(nDOF,nDOF);                   % stiffness matrix
Ktau  = zeros(nDOF,nDOF);                   % time stiffness matrix

F0    = zeros(nDOF,nDOFs);                  % time flux matrix for the initial condition
F1    = zeros(nDOF,nDOF);                   % time flux matrix at relative time tau=1

MDG    = zeros(nDOFs,nDOFs);                % global mass matrix
KxiDG  = zeros(nDOFs,nDOF);                 % global stiffness matrix
FRpDG  = zeros(nDOFs,nDOF);                 % flux on right interface on plus side
FRmDG  = zeros(nDOFs,nDOF);                 % flux on right interface on minus side
FLpDG  = zeros(nDOFs,nDOF);                 % flux on left interface on plus side
FLmDG  = zeros(nDOFs,nDOF);                 % flux on left interface on minus side


disp('Computing space-time integrals...')
%
% Compute space-time integrals
%
for itGP = 1:nGP
  for ixGP = 1:nGP
    xi  = xGPN(ixGP);
    tau = xGPN(itGP);
    weight = wGPN(itGP)*wGPN(ixGP);
    [theta,theta_xi,theta_tau]  = BaseFuncNodalTP_legendre(xi,tau);
    [psi,psi_xi]                = SpaceBaseFunc_2(xi);
    for k=[1:nDOF]
      for l=[1:nDOF]
        M(k,l)     = M(k,l)     + weight*theta(k)*theta(l);
        Kxi(k,l)   = Kxi(k,l)   + weight*theta(k)*theta_xi(l);
        Ktau(k,l)  = Ktau(k,l)  + weight*theta_tau(k)*theta(l);
      end
    end
    for k=[1:nDOFs]
      for l=[1:nDOF]
        KxiDG(k,l)  = KxiDG(k,l)  + weight*psi_xi(k)*theta(l);
      end
    end
  end
end
for ixGP = 1:nGP
  xi     = xGPN(ixGP);
  weight = wGPN(ixGP);
  [psi,psi_xi] = SpaceBaseFunc_2(xi);
  for k=[1:nDOFs]
    for l=[1:nDOFs]
      MDG(k,l) = MDG(k,l) + weight*psi(k)*psi(l);
    end
  end
end
iMDG = inv(MDG);


%
% Compute space integrals at t=1 and at t=0;
%
for ixGP = 1:nGP
    xi  = xGPN(ixGP);
    tau = 1;
    [theta1,theta1_xi,theta1_tau]=BaseFuncNodalTP_legendre(xi,tau);
    tau = 0;
    [theta0,theta0_xi,theta0_tau]=BaseFuncNodalTP_legendre(xi,tau);
    for k=[1:nDOF]
        for l=[1:nDOF]
            weight = wGPN(ixGP);
            F1(k,l)    = F1(k,l)    + weight*theta1(k)*theta1(l);
        end
    end
    [psi,psi_xi] = SpaceBaseFunc_2(xi);   % psi
    for k=[1:nDOF]
        for l=[1:nDOFs]
            weight = wGPN(ixGP);
            F0(k,l)    = F0(k,l)    + weight*theta0(k)*psi(l);
        end
    end
end

%
% Compute boundary time integrals
%
for itGP = 1:nGP
    weight = wGPN(itGP);
    tau = xGPN(itGP);
    xi  = 0;
    [theta0,theta0_xi,theta0_tau]   = BaseFuncNodalTP_legendre(xi,tau);
    [psi0,psi0_xi]                  = SpaceBaseFunc_2(xi);
    xi  = 1;
    [theta1,theta1_xi,theta1_tau]   = BaseFuncNodalTP_legendre(xi,tau);
    [psi1,psi1_xi]                  = SpaceBaseFunc_2(xi);
    for k=[1:nDOFs]
        for l=[1:nDOF]
            FRmDG(k,l)   = FRmDG(k,l)   + weight*psi1(k)*theta1(l);
            FRpDG(k,l)   = FRpDG(k,l)   + weight*psi1(k)*theta0(l);
            FLmDG(k,l)   = FLmDG(k,l)   + weight*psi0(k)*theta1(l);
            FLpDG(k,l)   = FLpDG(k,l)   + weight*psi0(k)*theta0(l);
        end
    end
end

%
K1  = F1 - Ktau;                    % The Balsara matrix
iK1 = inv(K1);                      % The inverse of the Balsara matrix

%%%%%%%%%%%%%%%%%%%%%%
% initial condition  %
%%%%%%%%%%%%%%%%%%%%%%
disp('Projecting initial condition...')
for i=[1:ghosts*2+IMAX]
    % map physical coordinates onto reference interval [0,1]
    for ixGP = 1:nGP
        j = j + 1;
        xi  = xGPN(ixGP);
        % periodic boundaries
        if (i <= ghosts)
            x_tmp = x(IMAX+i-1);
        elseif (i > ghosts+IMAX)
            x_tmp = x(i-IMAX-1);
        else
            x_tmp = x(i-ghosts);
        end
        % Recall xL=0; dx/2 is the determinant of the Jacobian matrix
        xGP = x_tmp - dx/2 + xi*dx;
        [psi,psi_xi]=SpaceBaseFunc_2(xi);
        %
        % (1) discontinuous initial condition
        % @Marten You want to try this
        % The initial condition is a Riemann problem
        % Ugly oscillations at the jump
        % Finite Volume could resolve this perfectly,
        % but DG s***s -- unless you use a limiter.
        if(xGP<0.5)
          u0 = 1;
        else
          u0 = 0;
        end
        %save the exact solution for later
        if (i > ghosts && i <= IMAX+ghosts)
          if(xGP<0.5)
            exact_solution(i-ghosts,ixGP) = 1;
          else
            exact_solution(i-ghosts,ixGP) = 0;
          end
        end


        % (2) Guass pulse, smooth
        %u0 = 1+exp(-(xGP-0.5)^2*20);
        %save the exact solution for later
        %if (i > ghosts && i <= IMAX+ghosts)
        %  exact_solution(i-ghosts,ixGP) = 1+exp(-(xGP-0.5)^2*20);
        %end


        % (3) periodic sin, smooth
        %u0 = 1+sin(2*pi()*xGP);
        %save the exact solution for later
        %if (i > ghosts && i <= IMAX+ghosts)
        %  exact_solution(i-ghosts,ixGP) = 1+sin(2*pi()*xGP);
        %end


        % (4) lone pulse, small discontinuity (undershoot)
        %u0 = max(0,exp(-(xGP-0.5).^2/.01)-10^-1);
        %save the exact solution for later
        % TODO
        % nicht sauber, da beim Plot die Unstetigkeit in der
        % analytischen Lsg nicht aufgeloest wird
        %if (i > ghosts && i <= IMAX+ghosts)
        %  exact_solution(i-ghosts,ixGP) = max(0,exp(-(xGP-0.5).^2/.01)-10^-1);
        %end


        % (5) ...
        % @Marten Feel free to try out other stuff. The advection equation
        % just transports you initial condition, so the exact_solution
        % is just a copy of your initial condition. Make sure your
        % example is periodic.


        % our initial DG polynomial
        uhat(:,i)  = uhat(:,i)  + wGPN(ixGP)*psi(:)*u0;
    end
    uhat(:,i) = iMDG*uhat(:,i);
end
%
disp('Preparation done. Starting code...')

%%%%%%%%%%%%%%%%%%%%%%
% main loop          %
%%%%%%%%%%%%%%%%%%%%%%
curTime  = 0;
while(curTime<tEnd)
  % Compute the timestep
  dt = CFL*CFLmax*dx/lmax;
  if(curTime+dt>tEnd)
      dt = tEnd-curTime;
  end

  % Local space-time discontinuous Galerkin predictor
  % Initial guess equal to old solution
  qhat = zeros(nDOF,IMAX+ghosts*2);  % predictor
  for i=[1:IMAX+ghosts*2]
    qhat(:,i) = uhat(1,i);
  end
  % Iterative local Space-Time DG method
  for iter=[1:N+2]
      for i=1:IMAX+ghosts*2
          % (1) Compute the flux f(Q)
          for iDOF = 1:nDOF
              fhat(iDOF,i) = (A*qhat(iDOF,i)')';
          end
          % (2) Compute the flux f*(Q) in the reference space
          fhatstar(:,i) = fhat(:,i)*dt/dx;
          % (3) Now do the ADER / local space-time DG iteration
          qhat(:,i) = iK1*( F0*uhat(:,i) - Kxi*fhatstar(:,i) );
      end
  end

  for i=ghosts+1:ghosts+IMAX
      %
      % Use a simple Rusanov flux at the element interfaces
      % (lmax is the maximum signal speed in the system)
      %
      FR = 0.5*FRmDG(:,:)*(fhat(:,i)+lmax*qhat(:,i)) + 0.5*FRpDG(:,:)*(fhat(:,i+1)-lmax*qhat(:,i+1));
      FL = 0.5*FLmDG(:,:)*(fhat(:,i-1)+lmax*qhat(:,i-1)) + 0.5*FLpDG(:,:)*(fhat(:,i)-lmax*qhat(:,i));
      %
      % The final arbitrary high order accurate DG scheme
      %
      uhat(:,i) = uhat(:,i) - iMDG*(dt/dx*(FR-FL) - KxiDG*fhatstar(:,i));
  end

  ## for i=2:IMAX-1
  ##      if( min(uhat(3,i)*uhat(3,i-1),uhat(3,i)*uhat(3,i+1))<0)
  ##          uhat(2,i) = minmod(uhat(1,i+1)-uhat(1,i),uhat(1,i)-uhat(1,i-1));
  ##          uhat(3:nDOFs,i) = 0;
  ##     end
  ## end

  %set periodic boundary
  for i=[1:ghosts]
    uhat(:,i) = uhat(:,i+IMAX+ghosts-1);
    uhat(:,i+ghosts+IMAX) = uhat(:,ghosts+i);
  end

  curTime = curTime + dt;


  % @Marten: if you want to speed up the simulation
  % do not plot every timestep but only e.g. every 0.1s
  % the timestep size is irregular and a scheme "plot
  % every 10th timestep" does not work here - I keep
  % it like this
  for i = 1:IMAX
    for ixGP = 1:nGP
      x_tmp = x(i);
      xi  = xGPN(ixGP);
      xGP = x_tmp - dx/2 + xi*dx;
      uplot(ixGP+(i-1)*nGP) = 0;
      [psi,psi_xi]=SpaceBaseFunc_2(xi);
      for k=1:nDOFs
          uplot(ixGP+(i-1)*nGP)  = uplot(ixGP+(i-1)*nGP)  + uhat(k,i+ghosts)*psi(k);
      end
      x_nodal(ixGP+(i-1)*nGP) = xGP;
    end
  end
  plot(x_nodal,uplot,'o');

  title(sprintf('t=%f',curTime))
  drawnow

end
%
% ==============================================================
% End of main loop
% ==============================================================
%
% Plot the exact solution
%
error = 0.0;
abs_val = 0.0;

for i = 1:IMAX
  for ixGP = 1:nGP
    x_tmp = x(i);
    xi  = xGPN(ixGP);
    xGP = x_tmp - dx/2 + xi*dx;
    u_solution(ixGP+(i-1)*nGP) = 0;
    [psi,psi_xi]=SpaceBaseFunc_2(xi);
    for k=1:nDOFs
      u_solution(ixGP+(i-1)*nGP)  = u_solution(ixGP+(i-1)*nGP)  + uhat(k,i+ghosts)*psi(k);
    end
  end
end

xCell = reshape(gaussPointsGlobal,N+1,cells);
yCell = reshape(u_solution,N+1,cells);

% TODO
% project nodel values onto equispaced
% grid to have a nice-looking plot
% or interpret as local Lagrange
% polynomials and plot

plot(x_nodal,u_solution,'o');
hold on;
plot(xCell,yCell,'linewidth',2);
hold on;
plot(x_nodal,reshape(exact_solution',IMAX*nGP,1),'r-');
title(sprintf('t=%f',curTime))
drawnow


% error computation
max_err = 0.0;
max_abs = 0.0;
error = 0.0;
abs_val = 0.0;

for i=[1:IMAX]
  for ixGP = 1:nGP
    cur_err = abs((exact_solution(i,ixGP)-u_solution(ixGP+(i-1)*nGP)));
    error = error + cur_err*cur_err;
    max_err = max(max_err,cur_err);
    cur_abs = abs(exact_solution(i,ixGP));
    abs_val = abs_val + cur_abs * cur_abs;
    max_abs = max(max_abs,cur_abs);
  end
end


disp('Linf Error was: ');
linf_error = max_err/max_abs;
disp(linf_error);

disp('Relative Error was: ');
rel_error = error/abs_val;
disp(rel_error);
