function [ Js, y, argmin, argmin2 ] = steady_functional_elliptic_boundary(z1, z2)

tic

%control domain
%\Gamma={0,1}.

%observation domain
%\omega_0=(0,1).

%weigthing parameter
%beta = 1/2.

w = [-10000:1:10000];
N = length(w);

%control-to-state map.
%cardinality of the partition of (0,1).

Nx = 48; 
% WARNING: Nx must be a multiple of 4, due to the definiton of the target z_d.
y = zeros(Nx,N);
Js = zeros(1,N);

% definition of discrete version of the target.
z_d = zeros(Nx,1);
for i = 1:(Nx/4)
    z_d(i) = z1;
    z_d((Nx/4)+i) = z2;
    z_d(2*(Nx/4)+i) = z2;
    z_d(3*(Nx/4)+i) = z1;
end

% computation of the control-to-state map.

% initialization i = 1.
zinit = zeros(Nx,1);
for i = 1:Nx
   %zinit(i,1) = (1-i/Nx)*w(1)+(i/Nx)*w(1);
   zinit(i,1) = w(1);
end
y(:,1) = poisson_cubic_nohomoboundary( Nx, zinit, w(1), w(1) );

for i = 2:N
    % we employ as initialization the state corresponding to the previous control.
    y(:,i) = poisson_cubic_nohomoboundary(Nx, y(:,i-1), w(i), w(i));
end

for i = 1:N
    % computation of the functional.
    Js(i) = (0.5)*(2*w(i).^2)+(0.5)^2*sum(((y(:,i)-z_d).^2)*(1/Nx));
end

% computation of the argmin.

Jsnotsorted = Js;

[Jssorted,indexes1] = sort(Js, 'ascend');

argmin = w(indexes1(1));
l = 2;
while (Jssorted(l-1) == Jssorted(l))
    argmin = [argmin;w(indexes1(l))];
    l = l+1;
end

% We compute
% argmin(2) = argmin(J_s \restriction_{[0,+\infty)}).

% SubStep 1 computation of the minimal index such that
% w(i)\geq 0.

i = 1;
while(w(i) < 0)
   i = i+1;
end

negativevaluescontrol = i-1;

Jsnegative = zeros(1,negativevaluescontrol);

for p = 1:negativevaluescontrol
   Jsnegative(p) = Jsnotsorted(p);
end

Jsnegativenotsorted = Jsnegative;

% SubStep 2 computation of argmin of
% argmin(J_s \restriction_{[0,+\infty)})

[Jsnegativesorted,indexes2] = sort(Jsnegative,'ascend');

argmin2 = w(indexes2(1));
l = 2;
while (Jsnegativesorted(l-1) == Jsnegativesorted(l))
    argmin2 = [argmin2;w(indexes2(l))];
    l = l+1;
end


level = (0.5)*(0.5)*(sum((z_d.^2)*(1/Nx)))*ones(N,1); % WARNING: Remember level = (\beta)*(0.5)*(sum((zdiscr.^2)*(1/Nx)))*ones(N,1);
figure(1);
clf(1);
%(30,144,255)
plot(w, Js, 'Color', 'r', 'LineWidth', 2.75);
%hold on;
%plot(w,level,'r','LineWidth',2);

%legend1 = legend("steady functional","$\frac12 \|z\|_{L^2}^2$",'Location','northeast');
%set(legend1,'Interpreter','latex');

%xlab1 = xlabel('u [control]','FontSize',10);
xlab1 = xlabel('control','FontSize', 10);
%set(xlab1,'Interpreter','latex');
%ylab1 = ylabel('$J_s$ [steady functional]','FontSize',10);
ylab1 = ylabel('functional','FontSize', 10);
%set(ylab1,'Interpreter','latex');
%xt = get(gca,'XTick');
set(gca,'FontSize',10);
ax = gca;
ax.LineWidth=1;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
exportgraphics(ax,'two_minimizers.pdf','ContentType','vector')

% figure(2);
% clf(2);
% plot(normxmin,Js);
% title('graph of the steady functional versus the norm of the state');
% xlabel('norm of y [state]','FontSize',20);
% ylabel('J_s [steady functional]','FontSize',20);
% xt = get(gca,'XTick');
% set(gca,'FontSize',20);
% figure(3);
% clf(3);
%     xone = xmin(1,:);
%     xtwo = xmin(2,:);
%     plot(xone,xtwo);
% title('representation of the set of controlled states [valid iff x_d\in \mathbb{R}^2]');
% xlabel('y_1 [first component]','FontSize',20);
% ylabel('y_2 [second component]','FontSize',20);
% xt = get(gca,'XTick');
% set(gca,'FontSize',20);

%set_plot_style_v03(36, 30, 6) % Prof. Sakamoto's way to set style of plots.

toc

end

