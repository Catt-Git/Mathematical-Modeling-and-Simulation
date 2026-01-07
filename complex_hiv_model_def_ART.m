function complex_hiv_model_def_ART
    clear; clc; close all;  
    
    % -----------------------------------
    % INITIAL CONDITION
    
    t_span = [0 5000];  % Days 0 to 14 years approx
    global t_start_art epsilon
    t_start_art = 300; % Days (Start therapy after acute phase)
    epsilon = 0.9;    % Efficacy (90% reduction in infection) 
    
    % LEGEND
    % 1: ST4n - Specific Heatlhy T4 non-effector
    % 2: ST4e - Specific Healthy  T4 effector
    % 3: ST4ni - Specific Latently infected T4 non-effector
    % 4: ST4ei - Specific Latently infected T4 effector
    % 5: T4n - Healthy Nonspecific T4 non-effector    
    % 6: T4e - Healthy Nonspecific T4 effector    
    % 7: T4ni - Nonspecific Latently infected T4 non-effector
    % 8: T4ei - Nonspecific Latently infected T4 effector
    % 9: ST8n - Healthy Specific T8 non-effector    
    % 10: ST8e - Healthy Specific T8 effector
    % 11: T8n - Healthy Nonspecific T8 non-effector
    % 12: T8e - Healthy Nonspecific T8 effector
    % 13: Tii - Actively infected T4 (viral burst)
    % 14: S - Antibodies
    % 15: Mi - Infected Macrophages
    % 16: V - Free Virus
    y0 = zeros(16,1);
    y0(1) = 0.92;   y0(2) = 0;      y0(3) = 0;      y0(4) = 0;   
    y0(5) = 920;    y0(6) = 379;    y0(7) = 0;      y0(8) = 0;   
    y0(9) = 0.47;   y0(10)= 0;      y0(11)= 467;    y0(12)= 150; 
    y0(13)= 0;      y0(14)= 0;      y0(15)= 0;      y0(16)= 3;  
    
    % The other parameters are below 
    % -----------------------------------------------------------------
    % SIMULATION
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'NonNegative', 1:16);
    [t, y] = ode15s(@ode_system, t_span, y0, options);
    
    % Process Results
    T4_total = sum(y(:,1:8), 2) + y(:,13); % Sum all T4 subsets
    T8_total = sum(y(:,9:12), 2); % Sum all T8 subsets
    V_RNA = y(:,16) * 2000; % Convert V to RNA
    V_RNA_scaled = V_RNA / 100; % Scale for plot
    Ratio = T8_total ./ (T4_total + eps);
    
% ----------------------------------------------
    % PLOT
    
    figure('Color', 'w', "Name", "Complex Model Dynamics with ART");
    
    % Figure 1
    subplot(2,2,1);
    plot_dynamics(t, T4_total, T8_total, V_RNA_scaled, [0 140]);
    title('Acute Infection Phase (0 - 140 Days)');
    
    % Figure 2
    subplot(2,2,2);
    plot_dynamics(t, T4_total, T8_total, V_RNA_scaled, [0 1000]);
    title('Chronic Phase (0 - 1000 Days)');
    
    % Figure 3
    subplot(2,2,3);
    plot_dynamics(t, T4_total, T8_total, V_RNA_scaled, [0 5000]);
    title('Progression to AIDS (0 - 5000 Days)');
    
    % Figure 4
    subplot(2,2,4);
    plot(t, Ratio, 'k-', 'LineWidth', 2);
    % Add ART line to Ratio plot
    if t(end) >= t_start_art
         xline(t_start_art, 'm--', 'Start ART');
    end
    xlabel('Days'); ylabel('Ratio');
    xlim([0 5000]); ylim([0 5.5]);
    grid on;
    title('T8/T4 Ratio');
end

function plot_dynamics(t, T4, T8, RNA, x_limits)
    global t_start_art
    plot(t, T4, 'b-', 'LineWidth', 1.5); hold on; 
    plot(t, T8, 'r-', 'LineWidth', 1.5);          
    plot(t, RNA, 'k--', 'LineWidth', 1.5);         
    
    % Add vertical line for ART start if within view
    if x_limits(2) >= t_start_art
        xline(t_start_art, 'm--', 'ART', 'LabelVerticalAlignment', 'bottom');
    end
    
    xlabel('Days');
    ylabel('Concentration'); % <--- Added Y-Axis Label Here
    xlim(x_limits);
    
    if x_limits(2) == 140
         ylim([0 2400]);
    else
         ylim([0 2200]);
    end
    grid on;
    legend('T4', 'T8', 'RNA/100', 'Location', 'Best');
end
    % --------------------------------------------------------
    % PARAMETERS
function dydt = ode_system(t, y)
    global t_start_art epsilon

    % Determine current drug efficacy
    if t >= t_start_art
        current_eps = epsilon;
    else
        current_eps = 0;
    end
       
    aM = 0.03;      % Conversion factor for the macrophage compartment 
    aT = 0.06;      % Conversion factor for the T cell compartment 
    c1 = 616.6;     % Michaelis-Menten half saturation for V in thymus
    c2 = 1000;      % Michaelis-Menten half saturation for V in blood
    delta = 0.001;  % Coefficient of reduced immune recognition of latently infected T4 and infected macrophages
    f = 0.524;      % T4/(T4+T8) ratio of neonate T cells
    k8 = 2.5;       % Rate T8e removes infected cells 
    km = 60;        % Rate macrophages remove virus 
    M_total = 360;  % Total macrophage density (healthy + infected) 
    mu_Ab = 0.023;  % Death rate of antibodies 
    mu_M = 0.1;     % Death rate of healthy and infected macrophages 
    mu_4e = 0.015;  % Death rate of T4 effector cells (T4e, ST4e, Ti4e, STi4e)
    mu_4n = 0.005;  % Death rate of T4 non-effector cells (T4n, ST4n, Ti4n, STi4n)
    mu_8e = 0.018;  % Death rate of T8 effector cells (T8e, ST8e)
    mu_8n = 0.006;  % Death rate of T8 non-effector cells (T8n, ST8n)
    mu_ii = 0.57;   % Death rate of actively infected cells (Tii)
    mu_V = 3.0;     % Death rate of free virus 
    N = 850;        % Number of virions per bursting Tii
    p = 0.03;       % Fraction/Probability of activated latently infected cells that become actively producing cells
    phi = 0.64;     % Fraction of T cells that differentiate
    pim = 32;       % Virus production rate by infected macrophages 
    rho_Ab = 1.55e5; % Maximal antibody production rate per helper T cell
    rho_4 = 1.98;   % Proliferation rate of ST4 in presence of virus 
    rho_8 = 0.36;   % Proliferation rate of ST8 in presence of infected cells
    r4 = 0.0097;    % Proliferation rate of T4 nonspecific cells 
    r8 = 0.0091;    % Proliferation rate of T8 nonspecific cells 
    nu = 0.001;     % ST4/Total T4 ratio in the thymic output
    % Time-varying parameters
    Age = 36;       % Age at inoculation 
    t1 = 7; t2 = 42; t3 = 40000; % Characteristic times for tropism/mutation
    
    % Helper Functions
    s_bar = 6.09;   % T cells flow from thymus to blood 
    s_t = s_bar * ((100 - Age - t/365)/100)^1.8;
    u_t = (t + t1) / (t + t2); % Simulates increasing tropism of virus for T cells
    
    % --- APPLY ART TO INFECTION TERMS ---
    kv_bar = 0.5;   % Baseline Rate of infection of T4 cells
    kv_t = kv_bar * (1 + t/t3) * u_t * (1 - current_eps); % Time-varying Virus infection rate
    
    e_bar = 0.064;  % Baseline Rate of infection of thymic T4 cells
    e_t = e_bar * (1 + t/t3) * u_t * (1 - current_eps); % Time-varying Thymic infection rate
    
    kvm_bar = 1.19; % Baseline Rate of infection of healthy macrophages
    kvm_t = kvm_bar * (1 + t/t3) * (1 - current_eps); % Time-varying Macrophage infection rate
    
    kAb_bar = 2.0e-7; % Baseline antibody efficacy/removal rate
    kAb_t = kAb_bar * (1 + t/t3) * u_t; % Time-varying antibody removal rate
    
    % State Variables
    ST4n=y(1); ST4e=y(2); ST4ni=y(3); ST4ei=y(4);
    T4n=y(5); T4e=y(6); T4ni=y(7); T4ei=y(8);
    ST8n=y(9); ST8e=y(10); T8n=y(11); T8e=y(12);
    Tii=y(13); S=y(14); Mi=y(15); V=y(16);
    
    % Algebraic Terms
    V_c1 = V / (c1 + V);
    V_c2 = V / (c2 + V);
    M_healthy = max(0, M_total - Mi); 
    Sum_Latent = ST4ni + T4ni + ST4ei + T4ei + (aT/aM)*Mi;
    R = rho_8 * (Tii + delta * Sum_Latent);
    
    % -------------------------------------------------------
    % ODE
    
    dydt = zeros(16,1);
    
    % T4 Specific
    dydt(1) = nu*f*s_t*(1 - e_t*V_c1) + (1-phi)*rho_4*V_c2*ST4n - kv_t*V_c2*ST4n - mu_4n*ST4n;
    dydt(2) = phi*rho_4*V_c2*ST4n - kv_t*V_c2*ST4e - mu_4e*ST4e;
    growth_factor = (1-phi)*(1-p) - p; 
    dydt(3) = nu*f*s_t*e_t*V_c1 + growth_factor*rho_4*V_c2*ST4ni + kv_t*V_c1*ST4n - (k8*delta*ST8e + mu_4n)*ST4ni;
    dydt(4) = phi*(1-p)*rho_4*V_c2*ST4ni + kv_t*V_c2*ST4e - (k8*delta*ST8e + mu_4e)*ST4ei;
    
    % T4 Nonspecific
    dydt(5) = (1-nu)*f*s_t*(1 - e_t*V_c1) + (1-phi)*r4*T4n - kv_t*V_c2*T4n - mu_4n*T4n;
    dydt(6) = phi*r4*T4n - kv_t*V_c2*T4e - mu_4e*T4e;
    dydt(7) = (1-nu)*f*s_t*e_t*V_c1 + growth_factor*r4*T4ni + kv_t*V_c2*T4n - (k8*delta*ST8e + mu_4n)*T4ni;
    dydt(8) = phi*(1-p)*r4*T4ni + kv_t*V_c2*T4e - (k8*delta*ST8e + mu_4e)*T4ei;
    
    % T8
    dydt(9) = nu*(1-f)*s_t + (1-phi)*R*ST8n - mu_8n*ST8n;
    dydt(10)= phi*R*ST8n - mu_8e*ST8e;
    dydt(11)= (1-nu)*(1-f)*s_t + (1-phi)*r8*T8n - mu_8n*T8n;
    dydt(12)= phi*r8*T8n - mu_8e*T8e;
    
    % Infected / Immune / Virus
    activation_source = p * (r4*T4ni + rho_4*V_c2*ST4ni);
    dydt(13)= activation_source - (k8*ST8e + mu_ii)*Tii;
    dydt(14)= rho_Ab*(ST4e + ST4ei)*V_c2 - kAb_t*V*S - mu_Ab*S;
    dydt(15)= kvm_t*V_c2*M_healthy - (k8*delta*ST8e + mu_M)*Mi;
    
    Prod_M = pim * Mi / aM;
    Prod_T = N * mu_ii * Tii / aT;
    Rem_M = km * V_c2 * (M_total / aM);
    Rem_Thymus = e_t * V_c1 * (s_t / aT);
    Total_T4_Healthy = ST4n + ST4e + T4n + T4e;
    Rem_Inf = (kv_t / aT) * V_c2 * Total_T4_Healthy;
    Rem_Ab = kAb_t * S * V;
    Rem_Nat = mu_V * V;
    
    dydt(16)= Prod_M + Prod_T - Rem_M - Rem_Thymus - Rem_Inf - Rem_Ab - Rem_Nat;
end