tic;

%%
n = 500;
dx = 1/(n-1); 

%%
Lap = diag(ones(n-1, 1),1) + diag(ones(n-1, 1),-1) - 2*eye(n);
Lap = (1/dx)^2*Lap;
    
% Heat
% A = Lap;

% Wave
A = zeros(2*n, 2*n); 
A(n+1:2*n,1:n) = Lap;
A(1:n,n+1:2*n) = eye(n);

%%
temp = zeros(1, n);
for i=1:n
    if (i>=n/4) && (i<=3*n/4)
        temp(1, i) = 1;
    end
end

% Heat
% B = diag(temp);

% Wave
% B = zeros(2*n,2*n);
% % B(n+1:2*n,n+1:2*n) = diag(temp);
% B(n+1:2*n,1:n) = diag(temp);

% Intro wave
B = zeros(2*n, 2*n);
B(n+1:2*n,1:n) = eye(n);
C = zeros(2*n, 2*n);
C(n+1:2*n,1:n) = Lap;


n_syst = size(A,1);
% rank(ctrb(A,B))
[Abar,Bbar,Cbar,T,k] = ctrbf(A,B,eye(n_syst)); 
nb_controllable_states = sum(k)

%% Hamiltonian matrix 

% Tracking
% Ham = [A, B*transpose(B); eye(n_syst), -transpose(A)];

% Intro wave
% Ham = [A, B; C, A];
% %Ham = [A, B; zeros(n_syst), A];

% No tracking
%Ham = [A, B*transpose(B); zeros(n_syst), -transpose(A)];

%% Spectrum
spectrum = eig(Ham);
re_spec = real(spectrum);
im_spec = imag(spectrum);
spec_abscissa = min(abs(real(spectrum)));

% Plots
clf()
plot(re_spec, im_spec, 'o','MarkerEdgeColor','k','MarkerFaceColor',[0.12 0.56 1], 'MarkerSize', 3);
hold on;
xh1 = xline(spec_abscissa, 'Color', [0.86 0.08 0.24], "LineWidth", 3, 'alpha', 1);
xh2 = xline(-spec_abscissa, 'Color', [0.86 0.08 0.24], "LineWidth", 3, 'alpha', 1);

% Set ConstantLine to 'back'
% Temporarily disable the warning that appears when accessing 
% the undocumented properties.
warnState = warning('off','MATLAB:structOnObject');
cleanupObj = onCleanup(@()warning(warnState)); 
Sxh1 = struct(xh1);         % Get undocumented properties (you'll get a warning)
Sxh2 = struct(xh2);         % Get undocumented properties (you'll get a warning)
clear('cleanupObj')         % Trigger warning reset
Sxh1.Edge.Layer = 'back';   % Set ConstantLine uistack
Sxh2.Edge.Layer = 'back';   % Set ConstantLine uistack

% Ticks
ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
grid minor
ax.GridLineStyle = ':';
ax.MinorGridLineStyle = ':';

% Saving 
% exportgraphics(ax,'eig_Ham_turnpike.pdf','ContentType','vector')
%exportgraphics(ax,'eig_Ham.pdf','ContentType','vector')
%exportgraphics(ax,'eig_Ham_intro.pdf','ContentType','vector')
exportgraphics(ax,'eig_Ham_turnpike_intro.pdf','ContentType','vector')

toc;