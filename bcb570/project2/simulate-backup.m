function simulate

% Parameters for Hes1 system
mu      = 3e-2;            % Degredation rate
Df      = 7.5e-10;         % Diffusion coefficient
alpha_p = 1.11e-2;         % Translation rate
tr_i    = 2e-12;           % Transcription initiation rate

% Spatial dimension
n   = 13;                  % Number of internal gridpoints
x   = linspace(0, 1, n+2)';% All gridpoints
x_i = x(2:end-1);          % Internal (non-boundary) gridpoints
h   = x(2) - x(1);         % Discrete spatial interval
n_m = 7;                   % Position of nuclear membrane (between grid points n_m-1 and n_m)

% Time dimension
dt  = 1e-2;                % Discrete time interval
t_f = 2;                   % Final time

% Construct basic A matrix
I = eye(n);                % Identity matrix
R = diag(ones(1,n-1),1);   % 1s on upper diagonal
A = R-2*I+R';              % Basic A matrix with -2s on diagonal and
                           % 1s on upper and lower diagonals

% Matrix for nuclear species
A_n = A;
A_n(n_m:end,:) = zeros(n-n_m+1,n);
A_n = A_n;
A_n(n_m-1,n_m) = 0;        % Conversion to cytoplasmic variant at nuclear membrane
M_n = I+dt*A_n/h^2;

% Matrix for cytoplasmic species
A_c = A;
A_c(1:n_m-1,:) = zeros(n_m-1,n);
A_c = A_c;
A_c(n_m,n_m-1) = 0;        % Conversion to cytoplasmic variant at nuclear membrane
A_c(end,end) = -1;         % Boundary condition
M_c = I+dt*A_c/h^2;

% Initialize vectors for nuclear and cytoplasmic variants of both species
r_n_0 = f(x);
r_n = r_n_0(2:end-1);
p_n_0 = f(x);
p_n = p_n_0(2:end-1);
r_c_0 = g(x);
r_c = r_c_0(2:end-1);
p_c_0 = g(x);
p_c = p_c_0(2:end-1);

% Run simulation
for tn = 1:ceil(t_f/dt)
  r_n = M_n*r_n + tr_i.*(1+(p_n/1e-2).^5);
  p_n = M_n*p_n;
  r_c = M_c*r_c;
  p_c = M_c*p_c + alpha_p*r_c;
  clf
  plot(x,r_c_0,'b:',x,r_n_0,'k:',x_i,r_c,'r.-',x_i,p_c,'g.-',x_i,r_n,'m.-',x_i,p_n,'y.-')
  title(sprintf('time t=%0.2f',tn*dt))
  drawnow
end

% Initial condition function for nuclear species
function y = f(x)
y = zeros(size(x));
q = x < 7;
y(q) = 1-(5*(x(q)+0.5).^2);

% Initial condition function for cytoplasmic species
function y = g(x)
y = zeros(size(x));
q = x >= 7;
y(q) = 1-(5*(x(q)-0.5).^2);
