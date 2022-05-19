clear
import casadi.*

%% Parameters:

% Space discretization for 1D
nx = 6;             
% Space discretization for 2D
NxSq = nx^2;         
% Time discretization
nt = 60;          
% Final time
T = 15;                 

%% Building a grid (Dirichlet problem!)
xline_ = linspace(0, 1, nx+2);
[xmsf,ymsf] = meshgrid(xline_, xline_);
xline = xline_(2:end-1);
[xms,yms] = meshgrid(xline, xline);
dx = xline(2) - xline(1);

%% Initial datum
Y0_aux = zeros(nx);
for i=1:nx
    for j=1:nx
        Y0_aux(i,j) = 0.1*sin(xline(i)*pi)*sin(xline(j)*pi);
    end
end
% Vectorizing 
Y0 = Y0_aux(:);

%% Dynamics matrices A and B

% B
B = zeros(nx^2, nx); 
B(1:nx, 1:nx) = eye(nx);        % control for i=1, ..., N_size.

% A
A1_mat = ones(nx, 1);
A2_mat = spdiags([A1_mat -2*A1_mat A1_mat], [-1 0 1], nx, nx);
I_mat  = speye(nx);
A_mat  = kron(I_mat, A2_mat) + kron(A2_mat, I_mat);
A = 0.02*A_mat;            

%% Time-stepping: CFL condition to change 'nt'.
tline = linspace(0, T, nt+1); 
dt = tline(2)-tline(1);

%% ---- Input variables ---------
opti = casadi.Opti();               %% CasADi function
X = opti.variable(NxSq, nt+1);      %% state trajectory
U = opti.variable(nx, nt+1);        %% control

%% ---- Dynamic constraints --------
M = eye(NxSq) - 0.5*dt/dx.^2*A;
L = eye(NxSq) + 0.5*dt/dx.^2*A;
P = 0.5*dt*B;

for k=1:nt 
   % Crank-Nicolson method
   opti.subject_to(M*X(:,k+1) == L*X(:,k) + 0.5*P*(U(:,k)+U(:,k+1)));
end

%% ---- State constraints --------
opti.subject_to(X(:,1)==Y0);
opti.subject_to(X(:,nt+1)==0);

%% ---- Optimization objective  ----------
Cost = (dx*sum(sum(U.^2))*(T/nt));
opti.minimize(Cost); 

%% ---- initial guesses for solver ---
Y = zeros(NxSq,nt+1);
opti.set_initial(X, Y);
opti.set_initial(U, 0);

%% ---- solve NLP ------
p_opts = struct('expand', true);
s_opts = struct('max_iter', 10000);     %% iteration limitation
opti.solver('ipopt', p_opts, s_opts);   %% set numerical backend

tic
sol = opti.solve();                     %% actual solve
Sol_u = sol.value(U);                   %% solved control function
Sol_x = sol.value(X);
time_axis = linspace(0, T, nt+1);

%% Plotting 
plot(time_axis, sum(Sol_u.^2)*dx, 'linewidth', 3, 'color', 'r')

ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
grid minor
exportgraphics(ax, 'heat_hum_control.pdf', 'ContentType', 'vector')
toc

% clf 
% Z = reshape(Sol_x(:,10), N_size, N_size);
% Z = [zeros(1, N_size+2) ; zeros(N_size, 1), Z, zeros(N_size, 1) ; zeros(1, N_size+2)];
% surf(xmsf, ymsf, Z)