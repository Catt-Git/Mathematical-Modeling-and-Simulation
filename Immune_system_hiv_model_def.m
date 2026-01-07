function Immune_system_hiv_model
    clear; clc; close all;  

    % -----------------------------------
    % PARAMETERS & CONFIGURATION
    t_final = 200; % Simulation duration (days) 
    p.t_start_art = 999;   % Day therapy starts
    p.epsilon_val = 0.90;  % Drug efficacy (90%)
    
    % Uninfected T Cells (T)
    p.lambda = 100;  % T-cell production rate        
    p.delta_T = 0.0046; % Natural T-cell death rate      
    
    % Infection Parameters
    p.beta = 4.8e-7; % Virus transmission 
    p.omega = 4.7e-7; % Cell-to-cell transmission 
    p.s = 0.01; % Saturation constant 
    
    % Infected T Cell Parameters (I)
    p.delta_I = 0.008;  % Infected cell death rate   
    p.m       = 0.001; % CTL neutralization rate     
    
    % Virus (V)
    p.N_vi    = 2000;   % Burst size (virions produced per infected cell)
    p.delta_V = 0.05;  % Viral clearance rate      
    p.b_z     = 0.01;      % Humoral neutralization rate 
    
    % Cytotoxic T Cell Response (E - CTL)
    p.p_I     = 0.002;  % CTL recruitment
    p.delta_E = 0.02;  % Death rate of CTL  
    
    
    % B Cell / Antibody Response (Z - Humoral)
    p.p_zv     = 0.0013;    % Humoral recruitment  
    p.delta_Z = 0.12;  % Humoral death rate     
    
    % --------------------------------------------------------
    % INITIAL CONDITIONS

    T0 = p.lambda / p.delta_T; 
    I0 = 1; 
    V0 = 0;
    E0 = 1;  
    Z0 = 1;  
    
    y0 = [T0; I0; V0; E0; Z0]; 
    
    % --------------------------------------------------------
    % SIMULATION

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
    [t, Y] = ode45(@(t,y) equations_lin_ode(t, y, p), [0 t_final], y0, opts);

    T = Y(:,1);      
    I = Y(:,2);  
    V = Y(:,3);      
    E = Y(:,4);      
    Z = Y(:,5);      
    
    v_copies_mL = V * 1000; 
    
    % ---------------------------------------------------------
    % PLOT

    figure('Name', 'Immune_HIV_model_ART', 'Color', 'w');
    
    % Subplot 1
    subplot(4, 1, 1);
    plot(t, T, 'b-', 'LineWidth', 2); hold on;
    % Only draw ART line if simulation went past start day
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'Start ART', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
    end
    ylabel({'Concentration'; '(cells/\muL)'}); 
    title('Uninfected CD4^+ T-Cells (T)');
    grid on; xlim([0 t_final]);
    
    % Subplot 2
    subplot(4, 1, 2);
    plot(t, I, 'r', 'LineWidth', 2); hold on;
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'LineWidth', 1.5);
    end
    ylabel({'Concentration'; '(cells/\muL)'}); 
    title('Infected CD4^+ T-Cells (I)');
    ylim([0 50]); 
    grid on; xlim([0 t_final]);
    
    % Subplot 3
    subplot(4, 1, 3);
    plot(t, E, '-', 'Color', [0.99, 0.05, 0.65], 'LineWidth', 2); hold on;
    plot(t, Z, 'Color', [0, 0.5, 1], 'LineWidth', 2);
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'LineWidth', 1.5);
    end
    ylabel({'Concentration'; '(cells/\muL)'}); 
    title('Immune Response (E & Z)');
    legend('CTL (E)', 'Antibodies (Z)', 'Location', 'Best');
    grid on; xlim([0 t_final]);
    
    % Subplot 4
    subplot(4, 1, 4);
    plot(t, v_copies_mL, 'k-', 'LineWidth', 2); hold on;
    if t_final >= p.t_start_art
        xline(p.t_start_art, 'm--', 'LineWidth', 1.5);
    end
    xlabel('Time (Days)');
    ylabel({'Viral Load'; '(copies/mL)'}); 
    title('Free Virus (V)');
    grid on; xlim([0 t_final]);
end

% ---------------------------------------------------------
% ODE SYSTEM

function dydt = equations_lin_ode(t, y, p)
    T = y(1);
    I = y(2); 
    V = y(3);
    E = y(4);
    Z = y(5);
    
    % Check if ART is active
    if t >= p.t_start_art
        current_epsilon = p.epsilon_val; % Drug active
    else
        current_epsilon = 0;             % Drug inactive
    end

    infection_rate = ((1 - current_epsilon) * p.beta * V * T) / (1 + p.s * V);
    

    dTdt = p.lambda - p.delta_T * T - infection_rate - p.omega * T * I;
    

    dIdt = infection_rate + p.omega * T * I - p.delta_I * I - p.m * E * I;
    

    dVdt = p.N_vi * p.delta_I * I - p.delta_V * V - p.b_z * V * Z;

    dEdt = p.p_I * I * E - p.delta_E * E;
    
    dZdt = p.p_zv * V * Z - p.delta_Z * Z;
    
    dydt = [dTdt; dIdt; dVdt; dEdt; dZdt];
end