clear all;
close all;
clc;
import casadi.*

alpha = 1;

%% Spatial discretization
li = 1;
Nx = 50;                    % number of interior points
hx = li/(Nx+1);             % space step
xaxe = linspace(0, li, Nx+2);

%% Initial datum
yi0 = sin(xaxe*pi);
yi1 = 2*exp(-(xaxe-li/2).^2)*0;

figure;
plot(xaxe, yi0, 'linewidth', 1.85, 'color', 'b')
title('Initial Datum');
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor

%% Setup 
ci = 1;
ht = hx/ci*0.9;             % time step
T = 5;                      % length of interval [0,T]
Nt = round(T/ht) - 1;       % number of interior points
taxe = linspace(0, T, Nt+2);

r = ht*ci/hx;               % for CFL: r <= 1

u = zeros(1, Nt+2);         % control = 0

%% Initialisation
yi = zeros(Nx+2, Nt+2);
yi(:, 1) = yi0;
yi(1, 2) =  ht*yi1(1) + (1-r^2)*yi0(1) + r^2*yi0(2);
yi(Nx+2, 2) = ht*yi1(Nx+2) + r^2*yi0(Nx+1) + (1-r^2)*yi0(Nx+2) + r^2*hx/ci*u(1);

for j=2:(Nx+1)
    yi(j, 2) = yi0(j) + ht*yi1(j) + 1/2*r^2*(yi0(j+1) - 2*yi0(j) + yi0(j-1));
end

%% Solve:
INx2 = ones(1, Nx +2); 
INx1 = ones(1, Nx +1);
A1 = r^2*diag(INx1, -1) + 2*(1-r^2)*diag(INx2) + r^2*diag(INx1, 1);
A1(1, 2) = 2*r^2;
A1(Nx+2, Nx+1) = 2*r^2;
B = zeros(Nx+2, 1);
B(Nx+2, 1) = 2*r^2*hx/ci;
A1 = sparse(A1);
B = sparse(B);

for k=3:(Nt+2)
    yi(:, k) = A1*yi(:, k-1) - yi(:, k-2) + B*u(k-1);
end

% %%%gif:
% figure;
% filename = 'WaveNeu_y_2ndOrder_beforeopti.gif';
% handle_line = plot(xaxe,yi0,'LineWidth',2); % plot ini datum
% for k=1:1:Nt+2
%     handle_line.YData = yi(:, k); % plot first part of z = (y, y_t)^T
%     axis([0,li,-1,6]);
%     xlabel('x'); ylabel('y(x)');
%     title('Wave equation: Neumann, y ODE45');
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

%% Casadi run        
opti = casadi.Opti();  % CasADi function
        
%% ---- Input variables ---------
Y = opti.variable(Nx+2,Nt+2); % state trajectory
U = opti.variable(1, Nt+2);   % control

%% ---- Control constraints -----------
% opti.subject_to(-100<=U(:));           % control is bounded
% opti.subject_to(U(:)<=100);            % control is bounded
opti.subject_to(Y(:,1)==yi(:,1));
%opti.subject_to(Y(:,2)==yi(:,2));
opti.subject_to(Y(:,Nt+2)==0);
opti.subject_to(Y(:,Nt+1)==0);
opti.set_initial(Y, yi);
opti.set_initial(U, 0);
 
%% ---- Dynamic constraints --------
opti.subject_to(Y(1, 2) ==  ht*yi1(1) + (1-r^2)*yi0(1) + r^2*yi0(2));
opti.subject_to(Y(Nx+2, 2) == ht*yi1(Nx+2) + r^2*yi0(Nx+1) + (1-r^2)*yi0(Nx+2) + r^2*hx/ci*U(1, 1));

for j=2:(Nx+1)
    opti.subject_to(Y(j, 2) == yi0(j) + ht*yi1(j) + 1/2*r^2*(yi0(j+1) - 2*yi0(j) + yi0(j-1)));
end

for k=3:(Nt+2)
    opti.subject_to(Y(:, k) == A1*Y(:, k-1) - Y(:, k-2) + B*U(1, k-1));
end

Ms = 2*eye(Nx+1); Ms(1, 1) = 1; Ms(Nx+1, Nx+1) = 1;
Mt = 2*eye(Nt+1); Mt(1, 1) = 1; Mt(Nt+1, Nt+1) = 1;
gint0 =  @(y) 1/hx^2*(y(2:Nx+2, Nt+1)-y(1:Nx+1, Nt+1)).^2*0; 
gint1 = @(y) alpha*Ms*gint0(y);
int1 = @(y, u) hx/2*sum( gint1(y) , 1) + u(1, 1:Nt+1).^2;
Func = @(y, u) ht/2*sum( int1(y, u)*Mt ) ;
opti.minimize(Func(Y, U));

%% ---- solve NLP              ------
opti.solver('ipopt'); % set numerical backend
tic
sol = opti.solve();   % actual solve
toc
 
state1 = opti.value(Y);
control1 = opti.value(U);
 
% %%%gif:
% figure;
% filename = 'WaveNeu_y_2ndOrder_afteropti.gif';
% handle_line = plot(xaxe,yi0,'LineWidth',2); % plot ini datum
% for k=1:1:Nt+2
%     handle_line.YData = state1(:, k); % plot first part of z = (y, y_t)^T
%     axis([0,li,-1,6]);
%     xlabel('x'); ylabel('y(x) opti');
%     title('Wave equation: Neumann, opti y ODE45');
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

figure;
plot(taxe(1:Nt+2), control1, 'linewidth', 1.85, 'color', 'b')
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
exportgraphics(ax,'wave_hum.pdf','ContentType','vector')
 
% F1=figure;

% plot(control1)
%         xticklabels({'0','5'});
%         xticks([0 Nt+1]);
%         xlabel('Time');
%         ylabel('Boundary Conrol');
%         line([0 2*Nt], [0 0],'Color','k');
%                 line([0 2*Nt], [1 1],'Color','k');
%                 xlim([0 Nt+2]);
%                               yticks([0 1])
%     yticklabels({'0','1'});
% set(gca,'FontSize',22)
%         a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',22)
%                set(F1,'PaperSize',[6 5]);
%         print(F1,'violation.pdf','-dpdf');
        
        

%             state=[flipud(state1); state1];
            
          