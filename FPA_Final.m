% Authors:  
%   Abdelrahman Sayed Sayed, <abdelrahman.ibrahim -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 25th of April 2025

%Neural Network Weight Matrices
A=[-1.20327 -0.07202 -0.93635; 1.18810 -1.50015 0.93519]; % `A` (2x3) 
B=[1.21464 -0.10502; 0.12023 0.19387; -1.36695 0.12201]; % `B` (3x2)
theta=[-7.5649*(10^-5); 1.34708*(10^-4); -6.24925*(10^-6)]; % Bias vector (3x1)
tau=10^6;
tau1=-(1/tau); %time constant for neurons
W=[[zeros(2),A]; [zeros(3,2),B*A]]; % composite weight matrix combining `A` and `B*A`

%Initial conditions for fixed point attractor and CTRNN
q0=[0.5;0.8]; % initial condition for the first two state (2x1)
s01=(B*q0)+theta; % remaining states (3x1)
s0=[q0; s01]; % Concatenates to form the full initial state vector (5x1)

% Initial set around this point
epsilon = 0.1*abs(s0);
x_low = s0-epsilon;
x_up = s0+epsilon;

% Definition of the function 
f = @(x,u) tau1*x+W*tanh(x); % $\dot{x} = f(x,u) = \tau x + W tanh(x)$
fpa = nonlinearSys('fpa',f);

n_x = 5; % state dimension
n_p = 1; % input dimension ~ not used in the code
p_low = 0; % input lower bound
p_up = 0; % input upper boud

%% Reach tube over approximation ~ Eq. to the 1st equation in proposition 2
% Parameters --------------------------------------------------------------

% Time range
t_init = 0;
t_final = 1;

params.tFinal = t_final;
params.R0 = zonotope(interval(x_low,x_up)); % x_in as initial state
params.U = interval(0,0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05; % Discretizes time into steps of 0.05 seconds
options.zonotopeOrder = 10; % Limits the complexity of zonotopes to avoid computational explosion
options.taylorTerms = 3; % Uses a third-order Taylor expansion for linearizing nonlinear dynamics
options.alg = 'lin';
options.tensorOrder = 2;
 
% Reachability Analysis ---------------------------------------------------
 
CORA_reach_output = reach(fpa, params, options); 

% Interval hull of all reach-tube steps
figure; hold on; % Fig.1
RT_interval_hull = interval(params.R0);
for i = 1:numel(CORA_reach_output.timeInterval.set)
    RT_interval_hull = RT_interval_hull | interval(CORA_reach_output.timeInterval.set{i}); % taking the union (|) with the interval of the current zonotope, forming the interval hull
    if i == 1
        % Plot the first zonotope and interval hull with legend labels
        plot(CORA_reach_output.timeInterval.set{i},[1,2],'b','DisplayName','Over-approximation of RT of all zonotopes') 
        plot(RT_interval_hull,[1,2],'r','DisplayName','Correction of each zonotope into interval') 
    else
        % Plot subsequent zonotopes and interval hulls without adding to legend
        plot(CORA_reach_output.timeInterval.set{i},[1,2],'b','HandleVisibility','off') 
        plot(RT_interval_hull,[1,2],'r','HandleVisibility','off') 
    end
end
RT_low = RT_interval_hull.inf; % reach tube lower bound
RT_up = RT_interval_hull.sup; % reach tube upper bound

% Compute successors from random initial states
sample_succ_number = 1000;
tic_sample_traj = tic;
rand_succ = NaN(n_x,sample_succ_number);

for i = 1:sample_succ_number
    x0 = x_low + rand(n_x,1).*(x_up-x_low);
    [~,x_traj] = ode45(@(t,x) f(x,0),[t_init t_final],x0); %using RK-4
    % Plot whole trajectories
    plot(x_traj(:,1),x_traj(:,2),'HandleVisibility','off') % Coloured Trajectories of 1000 sampled simulations of the nODE in Fig.1
    rand_succ(:,i) = x_traj(end,:)'; % These sampled final states are later used to visualize the system's output at t=1 (the black points in the last 2 figures)
end
toc_sample_traj = toc(tic_sample_traj)
set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$');
ylabel('$x_2$');
print('trajectories_and_zonotope_correction_fig.eps', '-depsc')

%% Error bounding by defining the function as a DT nonlinear system
%   and solving the reachability analylsis at tf=1s in a single step
%   (not intermediate time splitting)

% System Dynamics ---------------------------------------------------------


% Function f'*f/2
f_err = @(x,u) (tau1+W*diag(1-tanh(x).^2))*(tau1*x+W*tanh(x))/2;
error_function = nonlinearSysDT('error_function',f_err,t_final); 
 
% Parameters --------------------------------------------------------------
 
params.tFinal = t_final;
params.U = interval(0,0);
 
% Reachability Settings ---------------------------------------------------
clear options
options.zonotopeOrder = 10;
options.tensorOrder = 3;
options.errorOrder = 5;
 
% Reachability Analysis ---------------------------------------------------
% For each zonotope reach tube steps in the reachability analysis results
% we take the 1-step (no intermediate steps)reachable set over-approximation 
% for the discrete-time nonlinear system describing the static error function

figure % Fig.2 
hold on
grid on
OA_zonotopes = cell(numel(CORA_reach_output.timeInterval.set),1);
error_bound_DTreach = zeros(n_x,1);

for i = 1:numel(CORA_reach_output.timeInterval.set)
    % Define this intermediate zonotope reach-tube step as the initial set
    params.R0 = CORA_reach_output.timeInterval.set{i};
    % Do the 1-second DT reachability without intermediate steps
    R = reach(error_function, params, options);
    % Save the zonotope over-approximation output
    OA_zonotopes{i} = R.timePoint.set{2};
    % Plot this output OA, but exclude from legend
    plot(OA_zonotopes{i},[1,2],'HandleVisibility','off') 
    % Update error bound as error bound of union of each zonotope output
    error_bound_DTreach = max(error_bound_DTreach,max(abs(interval(OA_zonotopes{i})))); % (5x1)
end
error_bound_DTreach

% Interval hull of the error OA zonotopes
error_intHull = interval(OA_zonotopes{1});
for i = 1:numel(OA_zonotopes)
    error_intHull = convHull(error_intHull,interval(OA_zonotopes{i}));
end
plot(error_intHull,[1,2],'r','LineWidth',2,'DisplayName','$\Omega_{\varepsilon}$') 

set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
print('error_bound_fig.eps', '-depsc')

%% Main result on verification/reachability comparison
% Get nODE reach set based on ResNet output set and error bound
figure
hold on
grid on

input_set = taylm(interval(x_low,x_up));
plot(input_set+f(input_set),[1,2],'b','LineWidth',2,'DisplayName','$\Omega_{ResNet}$')
plot(input_set+f(input_set)+error_intHull,[1,2],'r','LineWidth',2,'DisplayName','$\Omega_{ResNet} + \Omega_{\varepsilon}$')

% Plot reachable set (at final time) of nODE from sample initial states
for j = 1:sample_succ_number
    plot(rand_succ(1,j),rand_succ(2,j),'k.','HandleVisibility','off') % random sampled points
end
% And its convex hull
nODE_reachset_convhull = convhull(rand_succ(1,:),rand_succ(2,:));
plot(rand_succ(1,nODE_reachset_convhull),rand_succ(2,nODE_reachset_convhull),'k','LineWidth',2,'DisplayName','$\mathcal{R}_{neural ODE}$')

% Safe set
saferegion = interval([0.2; 0.3], [0.6; 0.85]); 
plot(saferegion, [1, 2], 'g', 'LineWidth', 2.5, 'DisplayName','$\mathcal{X}_s$'); 

% legend('Location','northeast')
set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

% Set axis limits
xlim([0.15 0.65])
ylim([0.25 0.9])
print('ResNet_to_nODE_fig.eps', '-depsc')

%% Test in reverse
% Get resnet output set based on nODE reach set and error bound
figure
hold on
grid on

% Reach set of resnet belongs to output set of nODE - error bound
plot(CORA_reach_output.timePoint.set{end},[1,2],'k','LineWidth',2,'DisplayName','$\Omega_{neural ODE}$')

opposite_error = -error_intHull;
plot(CORA_reach_output.timePoint.set{end}+opposite_error,[1,2],'r','LineWidth',2,'DisplayName','$\Omega_{neural ODE} + \Omega_{-\varepsilon}$')

% Plot outputs of ResNet from random inputs
sample_succ_number = 1000;
rand_output = NaN(n_x,sample_succ_number);
for i = 1:sample_succ_number
    x0 = x_low + rand(n_x,1).*(x_up-x_low);
    rand_output(:,i) = x0+f(x0,0);
    plot(rand_output(1,i),rand_output(2,i),'b.','HandleVisibility','off')
end

% And its convex hull
resnet_output_convhull = convhull(rand_output(1,:),rand_output(2,:));
plot(rand_output(1,resnet_output_convhull),rand_output(2,resnet_output_convhull),'b','LineWidth',2,'DisplayName','$\mathcal{R}_{ResNet}$')

% Safe set
saferegion = interval([0.2; 0.3], [0.6; 0.85]); 
plot(saferegion, [1, 2], 'g', 'LineWidth', 2.5, 'DisplayName','$\mathcal{X}_s$'); 

set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

% Set axis limits
xlim([0.15 0.65])
ylim([0.25 0.9])
print('nODE_to_ResNet_fig.eps', '-depsc')