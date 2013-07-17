function simulate
% Daniel S. Standage
% May 1, 2012
% Adapted from code originally obtained from Benjamin Seibold at
% http://math.mit.edu/cse/codes/mit18086_fd_heateqn.m

% Parameters for Hes1 system
mu      = 3e-2;                % Degredation rate 
Df      = 7.5e-10;             % Diffusion coefficient
alpha_p = 1.11e-2;             % Translation rate
tr_i    = 2e-12;               % Transcription initiation rate

% Spatial dimension
n  = 13;                       % Number of internal grid points
x  = linspace(-1,1,n+2)';      % All grid points
xi = x(2:end-1);               % Internal (non-boundary) grid points
h  = x(2)-x(1);                % Discrete space interval
nm = ceil(n/2);                % First cytoplasmic grid point right of membrane

% Time dimension
dt = 1e-2;                     % Discrete time interval
tf = 20;                        % Final time

% Initialize vectors for nuclear and cytoplasmic variants of the mRNA and
% protein species
rn0 = g(x);                    % Nuclear mRNA
rn  = rn0(2:end-1);
pn0 = g(x);                    % Nuclear protein
pn  = pn0(2:end-1);   
rc0 = f(x);                    % Cytoplasmic mRNA
rc  = rc0(2:end-1);
pc0 = f(x);                    % Cytoplasmic protein
pc  = pc0(2:end-1);

I = eye(n);
R = diag(ones(1,n-1),1);

% M matrix for nuclear species
Dn = 2*I;
Dn(1,1) = 1;                   % Neumann boundary condition on left side
An =(R-Dn+R');
Mn = I-dt*mu+((Df*dt*An)/h^2); % Nuclear explicit time step
Mn(:,nm:end)=0;                % Nuclear membrane
Mn(nm:end,:)=0;

% M matrix for cytoplasmic species
Dc = 2*I;
Dc(end,end) = 1;                 % Neumann boundary condition on right side
Ac =(R-Dc+R');                              
Mc = I-dt*mu+((Df*dt*Ac)/h^2);   % Cytoplasmic explicit time step
Mc(:,1:nm)=0;                    % Nuclear membrane
Mc(1:nm,:)=0;

y_rn = zeros(1, ceil(tf/dt));
y_pn = zeros(1, ceil(tf/dt));
y_rc = zeros(1, ceil(tf/dt));
y_pc = zeros(1, ceil(tf/dt));

% Simulation
for tn = 1:ceil(tf/dt)
  %disp(sprintf('v=%dx%d, s=%dx%d\n', size(y_rn), size(rn)))
  rn = Mn*rn + tr_i.*(1+(pn/1e-2).^5); 
  pn = Mn*pn;
  pc = Mc*pc+alpha_p*rc;
  rc = Mc*rc;
  
  y_rn(tn) = sum(rn);
  y_pn(tn) = sum(pn);
  y_rc(tn) = sum(rc);
  y_pc(tn) = sum(pc);
  
  clf
  plot(x,rc0,'b:',x,rn0,'k:',xi,rc,'r.-',xi,pc,'g.-',xi,rn,'m.-',xi,pn,'y.-')
  title(sprintf('time t=%0.2f',tn*dt))
  drawnow
end

%clf
%plot(1:ceil(tf/dt), y_pc, 'b', 1:ceil(tf/dt), y_rc, 'r')

% Initial condition function for cytoplasmic species
function y = f(x)
y = zeros(size(x));
q = x > 0 ;
y(q) = 1-(5*(x(q)-0.5).^2);

% Initial condition function for nuclear species
function y = g(x)
y = zeros(size(x));
q = x < 0 ;
y(q) = 1-(5*(x(q)+0.5).^2);
