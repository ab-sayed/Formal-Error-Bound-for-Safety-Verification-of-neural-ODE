% Paper:
%   Bridging Neural ODE and ResNet: A Formal Error Bound for Safety Verification
%
% Authors:  
%   Abdelrahman Sayed Sayed, <abdelrahman.ibrahim -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%   Mohamed Ghazel, <mohamed.ghazel -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
%
% Date: 25th of April 2025
% Last update: 28th of April 2025
% Last revision: 28th of April 2025

%------------- BEGIN CODE --------------

figure
hold on
grid on

% Plot in the desired order: green, black, blue, red, pink
% Safe set (green)
saferegion = interval([0.2; 0.3], [0.6; 0.85]); 
plot(saferegion, [1, 2], 'g', 'LineWidth', 2.5, 'DisplayName','$\mathcal{X}_s$'); 

% Reachable set of Neural ODE (black)
for j = 1:sample_succ_number
    plot(rand_succ(1,j),rand_succ(2,j),'k.','HandleVisibility','off') % random sampled points
end
nODE_reachset_convhull = convhull(rand_succ(1,:),rand_succ(2,:));
plot(rand_succ(1,nODE_reachset_convhull),rand_succ(2,nODE_reachset_convhull),'k','LineWidth',2,'DisplayName','$\mathcal{R}_{neural ODE}$')

% ResNet output set (blue)
input_set = taylm(interval(x_low,x_up));
plot(input_set+f(input_set),[1,2],'b','LineWidth',2,'DisplayName','$\Omega_{ResNet}$')

% ResNet output set + your error bound (red)
plot(input_set+f(input_set)+error_intHull,[1,2],'r','LineWidth',2,'DisplayName','$\Omega_{ResNet} + \Omega_{\varepsilon}$')

% ResNet output set + Sander's error bound (pink)
plot(input_set+f(input_set)+error_sander,[1,2],'m','LineWidth',2,'DisplayName','$\Omega_{ResNet}$ + Error bound from Sander 2022')

set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

% Set axis limits
xlim([-0.5 1.6])
ylim([-0.5 1.6])
print('ResNet_to_nODE_comparison_with_Sander_fig.eps', '-depsc')


figure
hold on
grid on

% Plot in the desired order: green, blue, black, red, pink
% Safe set (green)
saferegion = interval([0.2; 0.3], [0.6; 0.85]); 
plot(saferegion, [1, 2], 'g', 'LineWidth', 2.5, 'DisplayName','$\mathcal{X}_s$'); 

% Reachable set of ResNet (blue)
sample_succ_number = 1000;
rand_output = NaN(n_x,sample_succ_number);
for i = 1:sample_succ_number
    x0 = x_low + rand(n_x,1).*(x_up-x_low);
    rand_output(:,i) = x0+f(x0,0);
    plot(rand_output(1,i),rand_output(2,i),'b.','HandleVisibility','off')
end
resnet_output_convhull = convhull(rand_output(1,:),rand_output(2,:));
plot(rand_output(1,resnet_output_convhull),rand_output(2,resnet_output_convhull),'b','LineWidth',2,'DisplayName','$\mathcal{R}_{ResNet}$')

% Neural ODE reach set (black)
plot(CORA_reach_output.timePoint.set{end},[1,2],'k','LineWidth',2,'DisplayName','$\Omega_{neural ODE}$')

% Neural ODE reach set + your error bound (red)
opposite_error = -error_intHull;
plot(CORA_reach_output.timePoint.set{end}+opposite_error,[1,2],'r','LineWidth',2,'DisplayName','$\Omega_{neural ODE} + \Omega_{-\varepsilon}$')

% Neural ODE reach set + Sander's error bound (pink)
plot(CORA_reach_output.timePoint.set{end}+error_sander,[1,2],'m','LineWidth',2,'DisplayName','$\Omega_{neural ODE}$ + Error bound from Sander 2022')

set(gca, 'FontSize', 12, 'FontName', 'Times')
lgd = legend('Location','best');
set(lgd, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

% Set axis limits
xlim([-0.5 1.6])
ylim([-0.5 1.6])
print('nODE_to_ResNet_comparison_with_Sander_fig.eps', '-depsc')

%------------- END OF CODE --------------