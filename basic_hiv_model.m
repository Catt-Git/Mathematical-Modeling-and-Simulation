function basic_hiv_model
    clear; clc; close all;
  
    % ---------------------------------------------------------------
    % PARAMETERS AND SETUP
    
    % Time parameters
    t_final = 800; % Total number of evaluated days 
    p.t_start_art = 200; % Day on which therapy would start
    
    % Model Parameters
    p.beta = 8e-7; % Infectivity rate
    p.NvI = 100; % Burst size (virions produced per infected cell)
    p.deltaV = 23; % Viral clearance rate
    p.lambda = 10000; % T-cell production rate 
    p.deltaT = 0.01; % Natural death rate of healthy T-cells
    p.deltaI = 0.7; % Death rate of infected T-cells
    
    % Therapy Efficacy
    p.epsilon_val = 0.9; % 90% efficacy
    
    % Initial Conditions 
    init_T_uL = 1000; % Initial T-cells per microliter
    init_I_uL = 0; % Initial Infected cells per microliter
    init_V = 5e4; % Initial Viral Load
    
    % Vector setup (conversion to mL for calculation)
    y0 = [init_T_uL * 1000; init_I_uL * 1000; init_V];
    
    % ----------------------------------------------------------
    %  SIMULATION
    
    p.rho = p.NvI * p.deltaI; % Calculate derived production rate 
    
    % Calc and print R0
    R0 = (p.beta * p.lambda * p.rho) / (p.deltaT * p.deltaI * p.deltaV);
    fprintf('Calculated R0 (Basic Reproduction Number): %.4f\n', R0);
    
    % Solver options
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    % Run ODE solver
    [t, y] = ode45(@(t,y) differential_equations(t, y, p), [0 t_final], y0, opts);
    
    % Prepare plot data (inverse conversion to uL for T and I)
    T_plot = y(:, 1) / 1000; 
    I_plot = y(:, 2) / 1000; 
    V_plot = y(:, 3);
    
    % =-----------------------------------------------------------
    % PLOT
    
    figure('Name', 'Basic_HIV_Model', 'Color', 'w');
    
    %Subplot 1:
    subplot(2, 1, 1);
    plot(t, T_plot, 'b-', 'LineWidth', 2); hold on;
    plot(t, I_plot, 'r--', 'LineWidth', 2);
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'Start ART', 'LabelVerticalAlignment', 'bottom');
        titleStr = 'CD4^+ Dynamics (With ART)';
    else
        titleStr = 'CD4^+ Dynamics (No ART)';
    end
    
    ylabel('Concentration (cells/\muL)');
    title(titleStr);
    legend('Uninfected T Cells (T)', 'Infected T Cells (I)', 'Location', 'Best'); 
    grid on;
    xlim([0 t_final]);
    
    % Subplot 2
    subplot(2, 1, 2);
    semilogy(t, V_plot, 'k-', 'LineWidth', 2); hold on;
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'Start ART');
    end
    
    ylabel('Viral Load (copies/mL)');
    xlabel('Time (Days)');
    title('Viral Load Dynamics');
    grid on;
    xlim([0 t_final]);
    
    % Set Y-axis limits to avoid errors if V drops to zero
    ylim([0.1 max(V_plot)*5]); 
end

    %-------------------------------------------------------------
    % ODE
function dydt = differential_equations(t, y, p)
    T = y(1); 
    I = y(2); 
    V = y(3);
    
    % Determine if therapy is active at current time 't'
    if t < p.t_start_art
        epsilon = 0;
    else
        epsilon = p.epsilon_val;
    end
    
    % Apply efficacy to parameter beta (infectivity)
    beta_eff = p.beta * (1 - epsilon);
    
    % Differential Equations
    dTdt = p.lambda - p.deltaT * T - beta_eff * V * T;
    dIdt = beta_eff * V * T - p.deltaI * I;
    dVdt = p.rho * I - p.deltaV * V;
    
    dydt = [dTdt; dIdt; dVdt]; 
end