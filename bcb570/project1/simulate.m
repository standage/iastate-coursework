% Daniel S. Standage
% BCB 570
% 26 Mar 2012
%
% Adapted from a script written by D.J. Higham, accessed at
% http://personal.strath.ac.uk/d.j.higham/chem/ssa_plot.m

clf
rand('state',100)

% Define stoichiometric matrix
V = [ 1   0  -1   0   0   0   0   0   0   0   0   0;
      0   1   0  -1   0   0   0   0   0   0   0   0;
      0   0   0   0   1   0  -1   0   0   0   0   0;
      0   0   0   0   0   1   0  -1   0   0   0   0;
      0   0   0   0   0   0   0   0   1   0  -1   0;
      0   0   0   0   0   0   0   0   0   1   0  -1; ];

% Constant parameters
tps_active = 0.5;    % strength of active promoter (transcripts / (promoter*sec) )
tps_repr   = 5.0E-4; % strength of repressed promoter (transcripts / (promoter*sec) )
tau_mRNA   = 2.0;    % mRNA half-life (min)
tau_prot   = 10.0;   % protein half-life (min)
K_m        = 40.0;   % (monomers/cell)
n          = 2;      % Hill coefficient

% Computed parameters
t_ave   = 2.89;   % average mRNA lifetime (min)
kd_mRNA = 0.347;  % mRNA decay rate (1/min)
kd_prot = 0.069;  % protein decay rate (1/min)
a_tr    = 29.97;  % active transcription rate (transcripts/min)
a0_tr   = 0.03;   % repressed transcription rate (transcripts/min)
k_tl    = 6.93;   % translation rate (proteins / (transcripts*min) )
alpha   = 216.4;  % (proteins / (promoters*cell*Km) )
alpha0  = 0.2164; % (proteins / (promoters*cell*Km) )
beta    = 0.2;    % ratio of mRNA to protein degredation rates

% Initial conditions for 6 species
X    = zeros(6,1);
X(1) = 20.0;

% Algorithm
t = 0;
tfinal = 1000;

count = 1;
tvals(1) = 0;
Xvals(:,1) = X;

while t < tfinal
     a(1) = a0_tr + a_tr*K_m^n / ( K_m^n+X(6)^n );
     a(2) = k_tl * X(1);
     a(3) = kd_mRNA * X(1);
     a(4) = kd_prot * X(2);
     a(5) = a0_tr + a_tr*K_m^n / ( K_m^n+X(2)^n );
     a(6) = k_tl * X(3);
     a(7) = kd_mRNA * X(3);
     a(8) = kd_prot * X(4);
     a(9) = a0_tr + a_tr*K_m^n / ( K_m^n+X(4)^n );
     a(10) = k_tl * X(5);
     a(11) = kd_mRNA * X(5);
     a(12) = kd_prot * X(6);

     asum = sum(a);
     j = min(find(rand<cumsum(a/asum)));
     tau = log(1/rand)/asum;
     X = X + V(:,j);

     count = count + 1;
     t = t + tau;
     tvals(count) = t;
     Xvals(:,count) = X;
end

%%%%%%%%%%% Plots

L = length(tvals);
tnew = zeros(1,2*(L-1));
tnew(1:2:end-1) = tvals(2:end);
tnew(2:2:end) = tvals(2:end);
tnew = [tvals(1),tnew];

Svals = Xvals(2,:);
ynew = zeros(1,2*L-1);
ynew(1:2:end) = Svals;
ynew(2:2:end-1) = Svals(1:end-1);
plot(tnew,ynew,'r-')
hold on

Svals = Xvals(4,:);
ynew = zeros(1,2*L-1);
ynew(1:2:end) = Svals;
ynew(2:2:end-1) = Svals(1:end-1);
plot(tnew,ynew,'Color',[1,0.7,0])

Svals = Xvals(6,:);
ynew = zeros(1,2*L-1);
ynew(1:2:end) = Svals;
ynew(2:2:end-1) = Svals(1:end-1);
plot(tnew,ynew,'b-')
hold off

legend('LacI protein','TetR protein','cI protein','Location','NorthWest');

xlabel('Time','FontSize',14)
ylabel('Species','FontSize',14)

axis([0 1000 0 6000])

set(gca,'FontWeight','Bold','FontSize',12)
grid on



