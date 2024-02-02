%% main file for solving optimal savings problem

clear
close all
clc;

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')

set(0,'DefaultTextFontSize', 12)
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultLineLineWidth',1)

temp = get(gca,'ColorOrder');
c1 = temp(1,:);
c2 = temp(2,:);
c3 = temp(3,:);

close all

%% define optimal savings structure

% preference and technology parameters
beta = 0.95; % discount factor
gamma = 1.5; % risk aversion
rf = 0.01; % risk-free rate
mu = 0.05; % expected return
sigma = 0.2; % volatility
theta = 0.5; % fraction of wealth invested in stocks
R = (1-theta)*exp(rf) + theta*exp(mu - sigma^2/2 + sigma*[1 -1]'); % gross return on wealth
Y = 1; % non-financial income
nz = length(R); % number of states
P = 0.5*ones(nz); % transition probability matrix

% define parameters for algorithm
N = 100; % number of grid points
N_plot = 1000; % number of grid points for plotting
MaxIter = 400; % maximum number of iterations
tol = 1e-5; % error tolerance

% grid for computation
aGrid = expGrid(0,100*Y,10*Y,N); % construct exponential grid

% grid for plotting
aMin = min(aGrid);
aMax = max(aGrid);
aGrid_plot = linspace(aMin,10*Y,N_plot);

%% model 1

os1.beta = beta*ones(nz);
os1.mu = @(c,z)(c^(-gamma)); % marginal utility
os1.R = repmat(R,1,nz);
os1.P = P; % transition probability matrix
os1.Y = Y*ones(nz);

rho = eigs(beta*os1.P.*os1.R,1);

os1.MaxIter = MaxIter;
os1.tol = tol;
os1.aGrid = aGrid;

% initialize consumption function
os1.Cmat0 = repmat(aGrid,nz,1);

tic
os1 = solve_os(os1);
toc

% plot results

% consumption function
figure
plot(os1.aGrid,os1.Cmat(1,:),'-','Color',c1);
xlabel('Asset')
ylabel('Consumption')
xlim([0 aMax])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os1_c','-dpdf')

% consumption function along iteration
imax = os1.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    c = interp1(os1.aGrid,os1.CMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,c,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Asset')
ylabel('Consumption')
xlim([0 10*Y])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os1_iter','-dpdf')

%% model 2

os2 = os1;
os2.mu = @(c,z)(exp(-gamma*c)); % marginal utility

tic
os2 = solve_os(os2);
toc

% plot results

% consumption function
figure
plot(os2.aGrid,os2.Cmat(1,:),'-','Color',c1);
xlabel('Asset')
ylabel('Consumption')
xlim([0 aMax])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os2_c','-dpdf')

% consumption function along iteration
imax = os2.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    c = interp1(os2.aGrid,os2.CMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,c,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Asset')
ylabel('Consumption')
xlim([0 10*Y])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os2_iter','-dpdf')

%% model 3

os3 = os1;

%C0 = Y*(3.5*(aGrid/(10*Y)).^2 + 0.01);
%C0 = min(C0,aGrid);

C0 = (sin(aGrid) + aGrid)/4;

%figure
%plot(aGrid,C0)
%xlim([0,10*Y])

os3.Cmat0 = repmat(C0,nz,1);

tic
os3 = solve_os(os3);
toc

% plot results

% consumption function
figure
plot(os3.aGrid,os3.Cmat(1,:),'-','Color',c1);
xlabel('Asset')
ylabel('Consumption')
xlim([0 aMax])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os3_c','-dpdf')

% consumption function along iteration
imax = os3.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    c = interp1(os3.aGrid,os3.CMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,c,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Asset')
ylabel('Consumption')
xlim([0 10*Y])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_os3_iter','-dpdf')
