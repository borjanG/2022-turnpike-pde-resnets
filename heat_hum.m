clear
%clc
import casadi.*

N_size = 6;     %% Space discretization for 1D
Nx = N_size^2;  %% Space discretization for 2D
Nt = 1000;       %% Time discretization
T = 15;         %% Final time

%% Discretization of the Space
xline_ = linspace(0,1,N_size+2);
[xmsf,ymsf] = meshgrid(xline_,xline_);
xline = xline_(2:end-1);
[xms,yms] = meshgrid(xline,xline);

dx = xline(2) - xline(1);

%% Initial data
Y0_2d = zeros(N_size);
for i=1:N_size
    for j=1:N_size
        Y0_2d(i,j) = 0.1*sin(xline(i)*pi)*sin(xline(j)*pi);
        %Y0_2d(i,j) = sin(xline(i)*pi);
    end
end
Y0 = Y0_2d(:);


%% Definition of the dynamics : Y' = AY+BU
B = zeros(N_size^2,N_size); 
B(1:N_size,1:N_size) = eye(N_size); %% control for i=1,...,N_size.

A1_mat = ones(N_size,1);
A2_mat = spdiags([A1_mat -2*A1_mat A1_mat], [-1 0 1], N_size,N_size);
I_mat  = speye(N_size);
A_mat  = kron(I_mat,A2_mat)+kron(A2_mat,I_mat);

A = 0.02*A_mat; %% A: Discrete Laplacian in 2D with uniform squared mesh

%% Discretization of the time : we need to check CFL condition to change 'Nt'.
tline = linspace(0,T,Nt+1); %%uniform time mesh

dt = tline(2)-tline(1);

%% Simulation of the uncontrolled trajectory
M = eye(Nx) - 0.5*dt/dx.^2*A;
L = eye(Nx) + 0.5*dt/dx.^2*A;
P = 0.5*dt*B;

Y = zeros(Nx,Nt+1);

Y(:,1) = Y0;
for k=1:Nt %% loop over time intervals
   %% Crank-Nicolson method without control
   Y(:,k+1) = M\L*Y(:,k);
end
YT = Y(:,Nt+1);

%% Target
ratio = 1.5;
%Y1 = ratio*YT; %% Target data
Y1_2d = zeros(N_size);
for i=1:N_size
    for j=1:N_size
        %Y0_2d(i,j) = sin(xline(i)*pi)*sin(xline(j)*pi);
        Y1_2d(i,j) = -2*xline(i)*exp(-xline(i)^2-xline(j)*pi);
    end
end
Y1 = Y1_2d(:);



%% 
% clf
% Z = reshape(Y(:,1),N_size,N_size);
% Z = [zeros(1,N_size+2) ; zeros(N_size,1), Z, zeros(N_size,1) ; zeros(1,N_size+2)];
% isurf = surf(xmsf,ymsf,Z,'FaceAlpha',0.3);
% isurf.CData = isurf.CData*0 + 10;
% hold on
% Z = reshape(Y(:,end),N_size,N_size);
% Z = [zeros(1,N_size+2) ; zeros(N_size,1), Z, zeros(N_size,1) ; zeros(1,N_size+2)];
% jsurf = surf(xmsf,ymsf,Z,'FaceAlpha',0.7);
% jsurf.CData = jsurf.CData*0 + 1;
% jsurf.Parent.Color = 'none';
% lightangle(10,10)
% legend({'Initial Condition','Final Data'})

opti = casadi.Opti();  %% CasADi function

%% ---- Input variables ---------
X = opti.variable(Nx,Nt+1); %% state trajectory
U = opti.variable(N_size,Nt+1);   %% control

%% ---- Dynamic constraints --------
for k=1:Nt %% loop over control intervals
   %% Crank-Nicolson method : this helps us to boost the optimization
   opti.subject_to(M*X(:,k+1)== L*X(:,k) + 0.5*P*(U(:,k)+U(:,k+1)));
end

%% ---- State constraints --------
opti.subject_to(X(:,1)==Y0);
%opti.subject_to(X(:,Nt+1)==0);

%% ---- Optimization objective  ----------
%Cost = (dx*sum(X(:,Nt+1)+3*ones(Nx,1))+dx*sum(sum(U.^2))*(T/Nt)+20*dx*sum(sum((X-2*ones(Nx,1)).^2))*(T/Nt));
Cost = (30*dx*sum(X(:,Nt+1)+10*ones(Nx,1))+dx*sum(sum(U.^2))*(T/Nt)+240*dx*sum(sum((X-Y1).^2))*(T/Nt));
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
Sol_u = sol.value(U); %% solved control function
Sol_x = sol.value(X);
time_axis = linspace(0, T, Nt+1);
%sum(Sol_u.^2)*dx
%sum(Sol_x.^2)*dx

plot(time_axis, sum(Sol_x.^2)*dx, 'linewidth', 3, 'color', 'b')
%plot(time_axis, sum(Sol_u.^2)*dx, 'linewidth', 3, 'color', 'r')

ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
%exportgraphics(ax,'heat_hum.pdf','ContentType','vector')

exportgraphics(ax,'heat_turnpike_y.pdf','ContentType','vector')
%exportgraphics(ax,'heat_turnpike_u.pdf','ContentType','vector')
toc


clf 
Z = reshape(Sol_x(:,10),N_size, N_size);
Z = [zeros(1,N_size+2) ; zeros(N_size,1), Z, zeros(N_size,1) ; zeros(1,N_size+2)];
surf(xmsf,ymsf,Z)
% Z = [zeros(1,N_size+2) ; zeros(N_size,1), Z, zeros(N_size,1) ; zeros(1,N_size+2)];
% isurf = surf(xmsf,ymsf,Z,'FaceAlpha',0.3);
% isurf.CData = isurf.CData*0 + 10;
% hold on
% Z = reshape(Sol_x(:,end),N_size,N_size);
% Z = [zeros(1,N_size+2) ; zeros(N_size,1), Z, zeros(N_size,1) ; zeros(1,N_size+2)];
% jsurf = surf(xmsf,ymsf,Z,'FaceAlpha',0.7);
% jsurf.CData = jsurf.CData*0 + 1;
% jsurf.Parent.Color = 'none';
% lightangle(10,10)
% legend({'Initial Condition','Final Data'})