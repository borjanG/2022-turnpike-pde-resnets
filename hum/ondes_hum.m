clear
clc
import casadi.*


%% We consider just one string i
% space
li = 1;
Nx = 100; % number of interior points
hx = li/(Nx+1); % space step
xaxe = linspace(0, li, Nx+2);

% time
T = 2; % length of interval [0,T]
Nt = 50; % number of interior points
ht = T/(Nt+1); % time step
taxe = linspace(0, T, Nt+2);

% data
ci = 20;
r = ht*ci/hx;
u = zeros(1, Nt+2);
%yi0 = sin(pi/2/li*xaxe);
%yi1 = pi/2/li*cos(pi/2/li*xaxe);
yi0 = 2*exp(-(xaxe-li/2).^2);
yi1 = 2*exp(-(xaxe-li/2).^2);

Y0 = [yi0, yi1];

%% Initialisation
yi = zeros(Nx+2, Nt+2);
yi(:, 1) = yi0;
yi(1, 2) =  ht*yi1(1) + (1-r^2)*yi0(1) + 1/2*r^2*yi0(2);
yi(Nx+2, 2) = ht*yi1(Nx+2) + 1/2*r^2*yi0(Nx+1) + (1-r^2)*yi0(Nx+2) + r^2*u(1);

for j=2:(Nx+1)
    yi(j, 2) = yi0(j) + ht*yi1(j) + 1/2*r^2*(yi0(j+1) - 2*yi0(j) + yi0(j-1));
end


%% Dirichlet Boudnary conditions
INx2 = ones(1, Nx +2); 
INx1 = ones(1, Nx +1);
A1 = r^2*diag(INx1, -1) + 2*(1-r^2)*diag(INx2) + r^2*diag(INx1, 1);
%A1(1, 2) = 2*r^2;
%A1(Nx+2, Nx+1) = 2*r^2;
%A2 = -eye(Nx+2);
B = zeros(Nx+2, 1);
B(1,1) = 1;

%% Simulation of the uncontrolled trajectory

Y = zeros(Nx,Nt+1);

opti = casadi.Opti();  %% CasADi function

%% ---- Input variables ---------
X = opti.variable(2*Nx,Nt+2); %% state trajectory
U = opti.variable(1, Nt+2);   %% control

%% ---- Dynamic constraints --------
for k=3:(Nt+2)
    opti.subject_to(X(:, k) == dt*A*X(:,k-1)-X(:,k-2)+B*u(k-1));
    yi(:, k) = A1*yi(:, k-1) - yi(:, k-2) + B*u(k-1);
end

%% ---- State constraints --------
opti.subject_to(X(:,1)==yi0);
opti.subject_to(X(:,2)==yi1);
opti.subject_to(X(:,Nt+1)==0);

%% ---- Optimization objective  ----------
Cost = (sum(U.^2)*(T/Nt));
opti.minimize(Cost); %% minimizing L2 over time

%% ---- initial guesses for solver ---
opti.set_initial(X, Y);
opti.set_initial(U, 0);

%% ---- solve NLP              ------
p_opts = struct('expand',true);
s_opts = struct('max_iter',10000); %% iteration limitation

opti.solver('ipopt',p_opts,s_opts); %% set numerical backend
tic
sol = opti.solve();   %% actual solve
Sol_x = sol.value(X);


Sol_u = sol.value(U); %% solved control function
time_axis = linspace(0, T, Nt+1);
plot(time_axis, Sol_u.^2, 'linewidth', 1.85, 'color', 'b')

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
exportgraphics(ax,'wave_hum.pdf','ContentType','vector')
toc