function cell_to_cell_hiv_model
    clear; clc; close all;
    
    % ----------------------------------------------------
    % PARAMETERS AND CONFIGURATION
    
    % Time and Therapy 
    t_final = 200; % Simulation duration (days) 
    p.enable_art  = false;        
    p.t_start_art = 200; % Day therapy starts
    p.epsilon = 0.90; % Drug efficacy (90%)
    
    % Biological Parameters
    p.lambda = 1e4; % T-cell production rate
    p.deltaT = 0.01; % Natural T-cell death rate
    p.beta = 2.4e-8; % Virus-to-Cell transmission rate
    p.omega = 1e-6; % Cell-to-Cell transmission rate 
    p.deltaI = 1.0; % Active infected cell death rate
    p.NvI = 2000; % Burst size (virions produced per infected cell)
    p.rho = p.NvI * p.deltaI; % Calculate derived production rate
    p.deltaV = 23.0; % Viral clearance rate
    p.a = 0.01; % Latency reactivation rate
    p.deltaL = 0.004; % Latent cell death rate
    p.alpha_frac = 0.001; % Fraction becoming latent from virus infection 
    p.eta = 0.001; % Fraction becoming latent from cell-to-cell 
    
    % --------------------------------------------------------
    % SIMULATION

    % Initial Conditions
    T0 = 1e6; % Initial Healthy T-cells (cells/mL)
    I0 = 0; % Initial Active Infected
    L0 = 0; % Initial Latent Infected
    V0 = 1e-3; % Initial Virus (small seed)
    
    y0 = [T0; I0; L0; V0];
    
    % R0 printing
    % This model distinguishes between free virus infection and direct cell contact
    R0_virus = (p.beta * T0 * p.rho) / (p.deltaI * p.deltaV);
    R0_cell  = (p.omega * T0) / p.deltaI;
    R0_total = R0_virus + R0_cell;
    fprintf('R0 Virus-to-Cell : %.4f\n', R0_virus);
    fprintf('R0 Cell-to-Cell  : %.4f\n', R0_cell);
    fprintf('R0 TOTAL         : %.4f\n', R0_total);

    if p.enable_art
        fprintf('Status           : ART ENABLED (Day %d)\n', p.t_start_art);
    else
        fprintf('Status           : ART DISABLED\n');
    end
    
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
    [t, y] = ode45(@(t,y) wang_model_equations(t, y, p), [0 t_final], y0, opts);
    
    % Unit Conversion (mL -> uL) for clinical plotting
    T = y(:, 1) / 1000;
    I = y(:, 2) / 1000;
    L = y(:, 3) / 1000;
    V = y(:, 4); 
    
    % ----------------------------------------------------------
    % PLOTTING
    
    figure('Name', 'CELL_TO_CELL_HIV_MODEL', 'Color', 'w');
    
    % Subplot 1
    subplot(3, 1, 1);
    plot(t, T, 'b-', 'LineWidth', 2); hold on;
    if p.enable_art && t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'Start ART', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
        titleStr = 'CD4^+ T-Cell Dynamics (With ART)';
    else
        titleStr = 'CD4^+ T-Cell Dynamics (No ART)';
    end
    
    ylabel('Concentration (cells/\muL)');
    title(titleStr);
    grid on; xlim([0 t_final]);
    
    % Subplot 2
    subplot(3, 1, 2);
    semilogy(t, I, 'r-', 'LineWidth', 1.5); hold on;
    semilogy(t, L, 'g-', 'LineWidth', 2);
    
    if p.enable_art && t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'LineWidth', 2);
    end
    
    ylabel('Infected Cells (cells/\muL)');
    legend('Active Infected (I)', 'Latent Infected (L)', 'Location', 'Best');
    title('Infected Cell Dynamics (Active vs Latent)');
    grid on; xlim([0 t_final]);
    ylim([1e-6 1000]); 
    
    % Subplot 3
    subplot(3, 1, 3);
    semilogy(t, V, 'k-', 'LineWidth', 2);
    
    if p.enable_art && t_final >= p.t_start_art
        xline(p.t_start_art, 'm--','LineWidth', 2);
    end
    
    ylabel('Viral Load (copies/mL)');
    xlabel('Time (Days)');
    title('Viral Load Dynamics');
    grid on; xlim([0 t_final]); 
    ylim([1e-2 max(V)*10]); 
end

    % ----------------------------------------------------
    % ODE
function dydt = wang_model_equations(t, y, p)
    T = y(1); 
    I = y(2); 
    L = y(3); 
    V = y(4);

    if p.enable_art && t >= p.t_start_art
        current_epsilon = p.epsilon;
    else
        current_epsilon = 0;
    end
    
    inf_virus = (1 - current_epsilon) * p.beta * V * T; 
    inf_cell  = p.omega * I * T;
    
    dTdt = p.lambda - p.deltaT * T - inf_virus - inf_cell;
    
    dIdt = (1 - p.alpha_frac) * inf_virus ...
         + (1 - p.eta)        * inf_cell ...
         + p.a * L ...
         - p.deltaI * I;
         
    dLdt = p.alpha_frac * inf_virus ...
         + p.eta        * inf_cell ...
         - p.a * L ...
         - p.deltaL * L;
         
    dVdt = p.rho * I - p.deltaV * V;
    
    dydt = [dTdt; dIdt; dLdt; dVdt];
end