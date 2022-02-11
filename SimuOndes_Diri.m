clear
clc

%% We consider just one string i
% space
li = 4*pi;
Nx = 100; % number of interior points
hx = li/(Nx+1); % space step
xaxe = linspace(0, li, Nx+2);

% time
T = 10; % length of interval [0,T]
Nt = 10000; % number of interior points
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

%figure;
%plot(xaxe, yi0);
%handle_line = plot(xaxe,yi0,'LineWidth',2);
%axis([0,li,-2,2]);
%xlabel('x'); ylabel('yi0');
%title('Wave equation');
%hold on;

%% Initialisation
yi = zeros(Nx+2, Nt+2);
yi(:, 1) = yi0;
yi(1, 2) =  ht*yi1(1) + (1-r^2)*yi0(1) + 1/2*r^2*yi0(2);
yi(Nx+2, 2) = ht*yi1(Nx+2) + 1/2*r^2*yi0(Nx+1) + (1-r^2)*yi0(Nx+2) + r^2*u(1);

for j=2:(Nx+1)
    yi(j, 2) = yi0(j) + ht*yi1(j) + 1/2*r^2*(yi0(j+1) - 2*yi0(j) + yi0(j-1));
end

%figure;
%plot(xaxe, yi(:, 1), xaxe, yi(:, 2));

%% Dirichlet Boudnary conditions
INx2 = ones(1, Nx +2); 
INx1 = ones(1, Nx +1);
A1 = r^2*diag(INx1, -1) + 2*(1-r^2)*diag(INx2) + r^2*diag(INx1, 1);
%A1(1, 2) = 2*r^2;
%A1(Nx+2, Nx+1) = 2*r^2;
%A2 = -eye(Nx+2);
B = zeros(Nx+2, 1);
for k=3:(Nt+2)
    yi(:, k) = A1*yi(:, k-1) - yi(:, k-2) + B*u(k-1);
end



% gif:
figure;
filename = 'SimuOndes_Diri.gif';
handle_line = plot(xaxe,yi0,'LineWidth',2);
for k=1:10:(Nt+2)
    handle_line.YData = yi(:,k);
    axis([0,li,-3,3]);
    xlabel('x'); ylabel('yi');
    title('Wave equation: Dirichlet');
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if k==1
           imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
