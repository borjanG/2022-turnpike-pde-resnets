clear
import casadi.*

Nx = 80;        %% Space discretization for 1D
T = 10;          %% Final time
Nt = T*Nx;      %% Time discretization

%% Discretization of the Space
xline_ = linspace(0,1,Nx+1);
xline = xline_(2:end);

dx = xline_(2) - xline_(1);

%% Initial data
Y0_1d = zeros(Nx, 1);
for i=1:Nx
    Y0_1d(i) = cos(xline(i)*pi);
end
Y0 = Y0_1d(:);


%% Definition of the dynamics : Y' = AY+BU
B = zeros(Nx,1); 
B(1) = 1/dx;

%Bvec = zeros(Nx,1);
%Bvec(25:Nx-1)=1;

% mol = zeros(Nx,1);
% for k=1:Nx
%     if xline(k)>0.25 && xline(k)<0.75
%         z = (2*xline(k)-0.25-0.75)/(0.5);
%         mol(k) = exp(-1/(1-z^2));
%     end
% end


% B = diag(mol);
%B=diag(Bvec);

A = -1*diag(ones(Nx,1))+diag(ones(Nx-1,1),-1);

%% Discretization of the time : we need to check CFL condition to change 'Nt'.
tline = linspace(0,T,Nt+1); %%uniform time mesh

dt = tline(2)-tline(1);
A = dt/dx*A;
B = dt*B;

opti = casadi.Opti();  %% CasADi function

%% ---- Input variables ---------
X = opti.variable(Nx,Nt+1); %% state trajectory
%U = opti.variable(Nx,Nt+1); %% control
U = opti.variable(Nt+1); %% control

%% ---- Dynamic constraints --------
for k=1:Nt %% loop over control intervals
   %opti.subject_to(X(:,k+1)== (eye(Nx)+A)*X(:,k) + B*(U(:,k)));
   opti.subject_to(X(:,k+1)== (eye(Nx)+A)*X(:,k) + B*(U(k)));
end

%% ---- State constraints --------
opti.subject_to(X(:,1)==Y0);
%opti.subject_to(U(2)-U(1)==-dt/dx*(Y0(2)-Y0(1)));
opti.subject_to(U(1)==Y0(1));
%opti.subject_to(X(:,Nt+1)==0);

%% ---- Optimization objective  ----------
%Cost = (30*dx*sum(X(:,Nt+1)-1.1*ones(Nx,1))+dx*sum(sum(U.^2))*(T/Nt)+240*dx*sum(sum((X-1).^2))*(T/Nt));
%Cost = (300*dx*sum(X(:,Nt+1)-1.5*ones(Nx,1))+sum(U.^2)*(T/Nt)+240*dx*sum(sum((X-1).^2))*(T/Nt));

Cost = (sum(U.^2)*(T/Nt)+1*dx*sum(sum((X-1).^2))*(T/Nt));

% M = eye(Nt+1);
% M(1,1)=0.5; M(Nt+1,Nt+1)=0.5;
% Z=X-0.5;
% M=T/Nt*M;
% Cost = (100000*sum(0.5*dx*(sum(Z(1:Nx-1,:).^2)+sum(Z(2:Nx,:).^2))*M)+(T/Nt)*(sum(U(1:Nt).^2)+sum(U(2:Nt+1).^2)));

opti.minimize(Cost); %% minimizing L2 over time

%% ---- initial guesses for solver ---
opti.set_initial(X, 0);
opti.set_initial(U, 0);

%% ---- solve NLP              ------
p_opts = struct('expand', true);
s_opts = struct('max_iter',10000); %% iteration limitation

opti.solver('ipopt',p_opts,s_opts); %% set numerical backend
tic
sol = opti.solve();   %% actual solve
Sol_u = sol.value(U); %% solved control function
Sol_x = sol.value(X);

% figure;
% filename = 'transport.gif';
% Yinit = zeros(Nx+1, 1);
% Yinit(1) = Sol_u(1);
% Yinit(2:Nx+1)=Y0;
% handle_line = plot(xline_, Yinit, 'LineWidth', 3.5); % plot ini datum
% Sol = zeros(Nx+1,Nt+1);
% Sol(2:Nx+1,:) = Sol_x;
% Sol(1,:) = Sol_u;
% ax = gca;
% ax.LineWidth=1.5;
% %ax.XGrid = 'on';
% %ax.YGrid = 'on';
% set(gca,'XMinorTick','on','YMinorTick','on')
% %grid minor
% for k=1:1:Nt+1
%     handle_line.YData = Sol(:, k);
%     axis([0,1,-1,6]);
%     xlabel('x'); ylabel('y_T(t,x)');
%     %title('Wave equation: Neumann, opti y ODE45');
%     drawnow;
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [AA,map] = rgb2ind(im,256);
%     if k==1
%         imwrite(AA,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
%     else
%         imwrite(AA,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%     end
% end

ax = gca;
%[tlist, xlist] = meshgrid(tline, xline);
[xlist, tlist] = meshgrid(xline_, tline);
%surf(tlist, xlist, Sol_x)
Sol = zeros(Nx+1,Nt+1);
Sol(2:Nx+1,:) = Sol_x;
Sol(1,:) = Sol_u;
surf(xlist, tlist, transpose(Sol))
xlabel('x')
ylabel('t')
view(2)
shading interp
colorbar
colormap(cool)
exportgraphics(ax,'transport_y.pdf','ContentType','vector')


% 
% time_axis = linspace(0, T, Nt+1);
% % 
% plot(time_axis, sum(Sol_x)*dx, 'linewidth', 3.5, 'color', 'b')
% %plot(time_axis, Sol_u, 'linewidth', 3.5, 'color', 'r')
% 
% ax = gca;
% ax.LineWidth=1.5;
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% set(gca,'XMinorTick','on','YMinorTick','on')
% %grid minor
% xlabel('time')
% %title('Control')
% title('State')
% % 
% exportgraphics(ax,'transport_turnpike_y.pdf','ContentType','vector')
% %exportgraphics(ax,'transport_turnpike_u.pdf','ContentType','vector')
% toc
