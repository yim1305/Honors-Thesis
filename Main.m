    %% Main
    %% By: Krittin Jindamai
clc; clear;

%% Orbital Parameters
mu = 3.986004418e14; 
Re = 6378137; % m         
a    = Re + 500e3;     
e = 0.01;
i = deg2rad(97.6);  % rad
RAAN = deg2rad(45); % rad
aop  = deg2rad(45); % rad
f    = 0;           % rad
n_orbit = sqrt(mu / a^3); % rad/s
orbit_param = [a,e,i,RAAN,aop,f];

%% Initial Orbit State
r0_PQW = (a*(1-e^2) / (1+(e*cos(f)))) * [cos(f);sin(f);0]; % m
r0_ECI = PQW2ECI(r0_PQW, orbit_param); % m
v0_PQW = (sqrt(mu/(a*(1-e^2)))) * [-sin(f);e+cos(f);0]; % m/s
v0_ECI = PQW2ECI(v0_PQW, orbit_param); % m/s
x_orbit = [r0_ECI; v0_ECI];

%% Simulation Time
T_one_orbit = 2*pi*sqrt(a^3/mu); % s
T = 15 * T_one_orbit;            % s
tspan = linspace(0,T,T*5);      
dt  = tspan(2) - tspan(1);

%% Initial Attitude State
q0   = eul2quat(0,0,0);        % [0;0;0;1]
w0_B = [0;0;-0];               % rad/s
x0_att = [q0; w0_B];    
x_att = x0_att;

%% Controller Parameters
[qd,wd_LVLH,wd_dot_LVLH,Rd] = compute_qd(r0_ECI,v0_ECI); % LVLH
R = quat2rotm_body(q0);
wd = R' * Rd * wd_LVLH; % body
wd_dot = R' * Rd * wd_dot_LVLH; % body

Kp = 0.3;
Kd = 1;

%% Actuator Parameters
Iw    = 42e-6;                 % kgm^2
theta = deg2rad(54.74);        % rad
delta_dot_max = deg2rad(60);   % deg/s 
Omega_dot_max = 55;            % rad/s2

% Initial gimbal angles
delta = [deg2rad(211); deg2rad(-31); deg2rad(211); deg2rad(-31)]; % rad

% flywheel speed parameter
Omega     = [838;838;838;838]; % rad/s (8000 RPM)
Omega_nom = [838;838;838;838]; % rad/s (8000 RPM)
Omega_max_scalar = 2000;       % rad/s
Omega_max = [Omega_max_scalar;Omega_max_scalar;Omega_max_scalar;Omega_max_scalar];
Omega_actual  = Omega;   

%% Steering Logic Parameters
lambda0 = 0.1;    mu1 = 5;
beta0   = 0.1;    mu2 = 5;    
alpha = 1;

%% Satellite Geometry  (12U / 230x240x360mm)
sat.Lx = 0.240;   % m 
sat.Ly = 0.230;   % m  
sat.Lz = 0.360;   % m 
sat.m  = 24;      % kg

% Center of Mass 
sat.CoM = [0; 0; 0.05];  % m

% Face definitions 
sat.faces(1).normal = [ 1;0;0]; sat.faces(1).area = sat.Ly*sat.Lz; sat.faces(1).cop = [ sat.Lx/2;0;0];
sat.faces(2).normal = [-1;0;0]; sat.faces(2).area = sat.Ly*sat.Lz; sat.faces(2).cop = [-sat.Lx/2;0;0];
sat.faces(3).normal = [ 0;1;0]; sat.faces(3).area = sat.Lx*sat.Lz; sat.faces(3).cop = [0; sat.Ly/2;0];
sat.faces(4).normal = [ 0;-1;0]; sat.faces(4).area = sat.Lx*sat.Lz; sat.faces(4).cop = [0;-sat.Ly/2;0];
sat.faces(5).normal = [ 0;0;1]; sat.faces(5).area = sat.Lx*sat.Ly; sat.faces(5).cop = [0;0; sat.Lz/2];
sat.faces(6).normal = [ 0;0;-1]; sat.faces(6).area = sat.Lx*sat.Ly; sat.faces(6).cop = [0;0;-sat.Lz/2];

% Drag coefficient 
sat.Cd = 2.2; 

% Inertia (J = diag([0.38, 0.37, 0.22]) ) kg/m2
Ixx_gc = sat.m/12 * (sat.Ly^2 + sat.Lz^2);
Iyy_gc = sat.m/12 * (sat.Lx^2 + sat.Lz^2);
Izz_gc = sat.m/12 * (sat.Lx^2 + sat.Ly^2);
J_gc   = diag([Ixx_gc, Iyy_gc, Izz_gc]);

d = sat.CoM;
Jc = J_gc - sat.m*(dot(d,d)*eye(3) - d*d'); % J about COM

% parallel axis for gimbals
m_cmg = 0.4; % kg
R_cmg1 = [0.01; 0; 0.1] - d;   % m
R_cmg2 = [-0.01; 0; 0.1] - d;  % m
R_cmg3 = [0; 0.01; 0.1] - d;   % m
R_cmg4 = [0; -0.01; 0.1] - d;  % m
R_cmg = [R_cmg1';R_cmg2';R_cmg3';R_cmg4'];
J = Jc;

for ii = 1:4
    r = R_cmg(ii, 1:3)'; 
    J_sum = m_cmg * ( (r' * r) * eye(3) - (r * r') );
    J = J + J_sum; 
end

%% MTQ parameter
sat.mtq.N         = 200;                           % turns
sat.mtq.I_max     = 0.4;                             % A
sat.mtq.A         = 0.02;                          % m2
sat.mtq.m_max     = sat.mtq.N * sat.mtq.I_max * sat.mtq.A;   % Am2 
sat.mtq.k_desat  = 3e-4;                                     % gain Am2/Nms 
sat.mtq.Omega_nom = Omega_nom;  % nominal flywheel speed rad/s

%% Scenario Switches
enable_mtq_desaturation = false; 

%% Storage Initialization
N = length(tspan);
x_att_store       = zeros(7, N);  
x_att_store(:,1)  = x0_att;
q_store           = zeros(4, N); 
q_store(:,1)      = q0;
wd_store          = zeros(3, N);  
wd_store(:,1)     = wd;
Omega_store       = zeros(4, N);  
Omega_store(:,1)  = Omega;
Omega_dot_store_actual  = zeros(4, N);
Omega_dot_store_cmd  = zeros(4, N);
delta_store       = zeros(4, N);  
delta_store(:,1)  = delta;
delta_dot_store_actual  = zeros(4, N);
delta_dot_store_cmd = zeros(4, N);
tau_cmd_store     = zeros(3, N);
tau_store         = zeros(3, N);
S_CMG_store       = zeros(1, N);
S_RW_store        = zeros(1, N);
alpha_store       = zeros(1, N);
Ee_store          = zeros(3, N);
ne_store          = zeros(1, N);
null_store        = zeros(4, N);
null_motion_store = zeros(4, N);
tau_dist_store    = zeros(3,N);
tau_gg_store      = zeros(3,N);
tau_aero_store    = zeros(3,N);
r_ECI_store       = zeros(3, N); 
r_ECI_store(:,1)  = r0_ECI;
x_orbit_store = zeros(6, N); 
x_orbit_store(:,1) = x_orbit;
tau_mtq_store    = zeros(3, N);
m_dipole_store   = zeros(3, N);
B_ECI_store      = zeros(3,N);
h_fw_store       = zeros(3,N);
h_excess_store   = zeros(3,N);
qd_store         = zeros(4,N);

%% Main Simulation Loop
for jj = 1:N-1

    %% 1. Extract current states
    t_current = tspan(jj);
    E_curr = x_att(1:3);
    n_curr = x_att(4);
    w_curr = x_att(5:7);
    q_curr = [E_curr; n_curr];
    R      = quat2rotm_body(q_curr);   % body -> ECI (current attitude)

    %% 2. Orbit state
    r_ECI     = x_orbit(1:3);
    v_ECI     = x_orbit(4:6);
    r_mag     = norm(r_ECI);
    r_ECI_mag = r_mag;
    nadir_I   = -r_ECI / r_mag;

    %% 3. Compute A, B, and S
    A     = matA(theta, delta);
    B     = matB(theta, delta);
    S_CMG = det(A*A');   
    S_RW  = det(B*B');

    %% 4. Angular momentum using Actual flywheel speed
    h_fw_B = Iw * B * Omega_actual;
    h_nom = Iw * B * Omega_nom; 
    h_max = Iw * B * Omega_max; 

    %% 5 Disturbances
    % Gravity gradient
    radial_I = r_ECI / r_ECI_mag;          
    o3_body  = ECI2BODY(q_curr, radial_I);
    tau_gg   = 3*mu/r_mag^3 * cross(o3_body, J*o3_body);
    
    % Aerodynamic drag
    [F_drag_body, a_drag, tau_aero] = compute_aero(q_curr, v_ECI, r_ECI, sat);
    
    % J2 
    a_J2 = compute_J2(x_orbit);
    
    % Total disturbance
    tau_dist = tau_gg + tau_aero;

    %% 6. Compute new qd, wd, wd_dot
    [qd,wd_LVLH,wd_dot_LVLH,Rd] = compute_qd(r_ECI,v_ECI); % LVLH

    % prevent qd from flipping sign midway
    if jj > 1
        if dot(qd, qd_prev) < 0
            qd = -qd;
        end
    end
    qd_prev = qd;

    wd     = R' * Rd * wd_LVLH;      
    wd_dot = R' * Rd * wd_dot_LVLH; 

    %% 7. Mode switch
    if t_current < 25
        alpha = 1;
    elseif t_current <= 35
        alpha = 0.5 * (1 + cos(pi*(t_current-25)/10));
    else
        alpha = 0;
    end

    %% 8. MTQ desaturation
    B_ECI = earth_magnetic_field2(x_orbit(1:3),t_current);
    B_body = R' * B_ECI;              % ECI to body

    h_excess = h_fw_B - h_nom;

    if enable_mtq_desaturation
        [tau_mtq, m_dipole] = compute_mtq(h_excess, sat, B_body);

        % MTQ always on while in RW mode (for now); off during CMG mode.
        if alpha ~= 0
            tau_mtq  = [0;0;0];
            m_dipole = [0;0;0];
        end
    else
        % Desaturation disabled (see enable_mtq_desaturation switch above).
        tau_mtq  = [0;0;0];
        m_dipole = [0;0;0];
    end

    %% 9. Eigen Axis Tracking Controller
    [tau_cmd, Ee, ne] = eigenaxis_controller(q_curr, qd, w_curr, wd_LVLH, wd_dot_LVLH, J, Kp, Kd, h_fw_B);

    %% 10. Steering law
    if alpha < 0.01
        d = zeros(4,1);  % no null motion in RW mode
    else
        d = null(A);
    end

    [lambda_val, beta_val, delta_dot_cmd, Omega_dot_cmd] = steering_logic(t_current, lambda0, mu1, S_CMG, ...
        beta0, mu2, S_RW, Iw, Omega_actual, alpha, A, B, tau_cmd, d);
 
    %% 11. Actuator saturation
    % Gimbal
    delta_dot_actual = max(-delta_dot_max, min(delta_dot_max, delta_dot_cmd));
    % Flywheel
    Omega_dot_actual = max(-Omega_dot_max, min(Omega_dot_max, Omega_dot_cmd));

    %% 12. Actual torque from actuator outputs
    tau = (Iw * A * diag(Omega_actual) * delta_dot_actual) + (Iw * B * Omega_dot_actual);
    tau_total = tau - tau_mtq - tau_dist;
        
    %% 13. Attitude dynamics with actual torque
    x_att_dot = attitude_dynamics(x_att, J, tau_total, h_fw_B);

    %% 14. Integrate attitude (quaternion)
    E_next = E_curr + x_att_dot(1:3) * dt;
    n_next = n_curr + x_att_dot(4)   * dt;
    w_next = w_curr + x_att_dot(5:7) * dt;

    % Renormalize quaternion
    q_norm = norm([E_next; n_next]);
    E_next = E_next / q_norm;
    n_next = n_next / q_norm;
    x_att  = [E_next; n_next; w_next];

    %% 15. Update gimbal angles and flywheel speed
    delta  = delta + delta_dot_actual * dt; 
    Omega_actual = Omega_actual + Omega_dot_actual * dt;
    Omega_actual = max(-Omega_max, min(Omega_max, Omega_actual));  
    Omega = Omega_actual;

    %% 16. Update Orbit
    x_orbit_dot = orbit_dynamics(x_orbit, a_drag, a_J2);
    x_orbit(4:6) = x_orbit(4:6) + x_orbit_dot(4:6) * dt;
    x_orbit(1:3) = x_orbit(1:3) + x_orbit(4:6) * dt;

    %% 17. Store at jj+1
    q_store(:,jj+1)          = [E_next; n_next];
    x_att_store(:,jj+1)      = x_att;
    wd_store(:,jj+1)         = wd;
    Omega_store(:,jj+1)      = Omega;
    Omega_dot_store_actual(:,jj+1)  = Omega_dot_actual;
    Omega_dot_store_cmd(:,jj+1)  = Omega_dot_cmd;
    delta_store(:,jj+1)      = delta;
    delta_dot_store_actual(:,jj+1)  = delta_dot_actual;
    delta_dot_store_cmd(:,jj+1)  = delta_dot_cmd;
    tau_cmd_store(:,jj+1)    = tau_cmd;
    tau_store(:,jj+1)        = tau;
    S_CMG_store(jj+1)        = S_CMG;
    S_RW_store(jj+1)         = S_RW;
    alpha_store(jj+1)        = alpha;
    Ee_store(:,jj+1)         = Ee;
    ne_store(jj+1)           = ne; 
    null_store(:,jj+1)       = d;
    null_motion_store(:,jj+1) = beta_val * d;
    tau_dist_store(:,jj+1) = tau_dist;
    tau_gg_store(:,jj+1) = tau_gg;
    tau_aero_store(:,jj+1) = tau_aero;
    x_orbit_store(:,jj+1) = x_orbit;
    tau_mtq_store(:,jj+1)  = tau_mtq;
    m_dipole_store(:,jj+1) = m_dipole;
    B_ECI_store(:,jj+1) = B_ECI;
    h_excess_store(:,jj+1) = h_excess;
    h_fw_store(:,jj+1)     = h_fw_B;
    qd_store(:,jj+1) = qd;
end

%% Plots
lw = 1.5; % consistent line width across every figure

figure(1); % Attitude
tiledlayout(2,1,'TileSpacing','compact');

nexttile; % Error quaternion
plot(tspan, Ee_store(1,:), 'r-', tspan, Ee_store(2,:), 'g-', ...
     tspan, Ee_store(3,:), 'b-', tspan, ne_store, 'k-', 'LineWidth', lw);
grid on; ylabel('Error Quaternion Component [-]');
title('Attitude Error Quaternion');
legend('\epsilon_{e1}','\epsilon_{e2}','\epsilon_{e3}','\eta_e');
ylim([-1.1 1.1]);

nexttile; % Euler angles
euler_store = rad2deg(quat2eul(q_store));
plot(tspan, euler_store(1,:), 'r-', tspan, euler_store(2,:), 'g-', ...
     tspan, euler_store(3,:), 'b-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Euler Angle [deg]');
title('Attitude Euler Angles');
legend('roll','pitch','yaw');

figure(2); % Angular Velocity
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(tspan, rad2deg(x_att_store(5,:)), 'r-', ...
     tspan, rad2deg(x_att_store(6,:)), 'g-', ...
     tspan, rad2deg(x_att_store(7,:)), 'b-', 'LineWidth', lw);
grid on; ylabel('Angular Velocity [deg/s]');
title('Satellite Angular Velocity');
legend('\omega_x','\omega_y','\omega_z');

nexttile;
plot(tspan, rad2deg(x_att_store(5,:)-wd_store(1,:)), 'r-', ...
     tspan, rad2deg(x_att_store(6,:)-wd_store(2,:)), 'g-', ...
     tspan, rad2deg(x_att_store(7,:)-wd_store(3,:)), 'b-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Angular Velocity Error [deg/s]');
title('Satellite Angular Velocity Error');
legend('\omega_x','\omega_y','\omega_z');

figure(3); % Gimbal
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(tspan, rad2deg(delta_dot_store_actual(1,:)), 'r-', ...
     tspan, rad2deg(delta_dot_store_actual(2,:)), 'g-', ...
     tspan, rad2deg(delta_dot_store_actual(3,:)), 'b-', ...
     tspan, rad2deg(delta_dot_store_actual(4,:)), 'k-', 'LineWidth', lw);
grid on; ylabel('Gimbal Rate [deg/s]');
title('Gimbal Rates (Actual)');
legend('\delta_1','\delta_2','\delta_3','\delta_4');

nexttile;
plot(tspan, rad2deg(delta_store(1,:)), 'r-', ...
     tspan, rad2deg(delta_store(2,:)), 'g-', ...
     tspan, rad2deg(delta_store(3,:)), 'b-', ...
     tspan, rad2deg(delta_store(4,:)), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Gimbal Angle [deg]');
title('Gimbal Angles');
legend('\delta_1','\delta_2','\delta_3','\delta_4');

figure(4); % Reaction Wheel
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(tspan, Omega_dot_store_actual(1,:), 'r-', ...
     tspan, Omega_dot_store_actual(2,:), 'g-', ...
     tspan, Omega_dot_store_actual(3,:), 'b-', ...
     tspan, Omega_dot_store_actual(4,:), 'k-', 'LineWidth', lw);
grid on; ylabel('Flywheel Acceleration [rad/s^2]');
title('Flywheel Acceleration');
legend('\Omega_1','\Omega_2','\Omega_3','\Omega_4');

nexttile;
plot(tspan, Omega_store(1,:), 'r-', ...
     tspan, Omega_store(2,:), 'g-', ...
     tspan, Omega_store(3,:), 'b-', ...
     tspan, Omega_store(4,:), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Flywheel Speed [rad/s]');
title('Flywheel Speed');
legend('\Omega_1','\Omega_2','\Omega_3','\Omega_4');

figure(5); % Torque HCMG
plot(tspan, tau_cmd_store(1,:)*1000, 'r--', ...
     tspan, tau_cmd_store(2,:)*1000, 'g--', ...
     tspan, tau_cmd_store(3,:)*1000, 'b--', ...
     tspan, tau_store(1,:)*1000,     'r-',  ...
     tspan, tau_store(2,:)*1000,     'g-',  ...
     tspan, tau_store(3,:)*1000,     'b-',  'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Torque [mNm]');
title('Commanded vs Actual HCMG Torque');
legend('\tau_{cmd,x}','\tau_{cmd,y}','\tau_{cmd,z}', ...
       '\tau_{act,x}','\tau_{act,y}','\tau_{act,z}');

figure(6); % Sing Params
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(tspan, S_CMG_store, 'b-', tspan, S_RW_store, 'r-', 'LineWidth', lw);
grid on; ylabel('Singularity Parameter [-]');
title('Singularity Parameters');
legend('S_{CMG}','S_{RW}');

nexttile; % Null motion
plot(tspan, rad2deg(null_motion_store(1,:)), 'r-', ...
     tspan, rad2deg(null_motion_store(2,:)), 'g-', ...
     tspan, rad2deg(null_motion_store(3,:)), 'b-', ...
     tspan, rad2deg(null_motion_store(4,:)), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Null Motion Rate [deg/s]');
title('Null Motion');
legend('HCMG 1','HCMG 2','HCMG 3','HCMG 4');

figure(7); % Mode Switch Param
plot(tspan, alpha_store, 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('\alpha [-]');
title('Mode Switch Parameter \alpha');
ylim([-0.1 1.1]);

figure(8); % Disturbances
tiledlayout(3,1,'TileSpacing','compact');

nexttile;
plot(tspan, tau_dist_store(1,:)*1e6,'r-', tspan, tau_dist_store(2,:)*1e6,'g-', ...
     tspan, tau_dist_store(3,:)*1e6,'b-','LineWidth',lw);
grid on; ylabel('Torque [\muNm]'); title('Total Disturbance');
legend('\tau_x','\tau_y','\tau_z');

nexttile;
plot(tspan, tau_gg_store(1,:)*1e6,'r-', tspan, tau_gg_store(2,:)*1e6,'g-', ...
     tspan, tau_gg_store(3,:)*1e6,'b-','LineWidth',lw);
grid on; ylabel('Torque [\muNm]'); title('Gravity Gradient');
legend('\tau_x','\tau_y','\tau_z');

nexttile;
plot(tspan, tau_aero_store(1,:)*1e6,'r-', tspan, tau_aero_store(2,:)*1e6,'g-', ...
     tspan, tau_aero_store(3,:)*1e6,'b-','LineWidth',lw);
grid on; xlabel('Time [s]'); ylabel('Torque [\muNm]'); title('Aerodynamic Drag');
legend('\tau_x','\tau_y','\tau_z');

figure(9); % flywheel angular moomentum + mtq
h_excess_mag_store = vecnorm(h_excess_store);

subplot(4,1,1); % excess h
plot(tspan, h_excess_store(1,:), 'r-', ...
     tspan, h_excess_store(2,:), 'g-', ...
     tspan, h_excess_store(3,:), 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('[kg\cdotm^2/s]');
title('Excess Angular Momentum');
legend('h_x','h_y','h_z');

subplot(4,1,2)
plot(tspan,h_excess_mag_store,'k','LineWidth',2);
title('Excess Angular Momentum Magnitude');
grid on; xlabel('Time [s]'); ylabel('[kg\cdotm^2/s]');
hold on
%yline(0.7*norm(h_max),'g--','70%');
yline(norm(h_max),'k--','h_{max}');

subplot(4,1,3); % dipole moment
plot(tspan, m_dipole_store(1,:), 'r-', ...
     tspan, m_dipole_store(2,:), 'g-', ...
     tspan, m_dipole_store(3,:), 'b-', 'LineWidth', 1.5);
yline(sat.mtq.m_max,  'k--', 'Max');
yline(-sat.mtq.m_max, 'k--');
grid on; xlabel('Time [s]'); ylabel('[Am²]');
title('MTQ Dipole Moment Command');
legend('m_x','m_y','m_z');

subplot(4,1,4); % desat torque
plot(tspan, tau_mtq_store(1,:)*1e6, 'r-', ...
     tspan, tau_mtq_store(2,:)*1e6, 'g-', ...
     tspan, tau_mtq_store(3,:)*1e6, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('[\muNm]');
title('MTQ Desaturation Torque');
legend('\tau_{mtq,x}','\tau_{mtq,y}','\tau_{mtq,z}');

%% 3D
figure(10);
plot3(x_orbit_store(1,:)/1000, x_orbit_store(2,:)/1000, x_orbit_store(3,:)/1000, ...
      'b-', 'LineWidth', 2);
hold on;
grid on;
plot3(x_orbit_store(1,1)/1000, x_orbit_store(2,1)/1000, x_orbit_store(3,1)/1000, ...
      'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
plot3(x_orbit_store(1,end)/1000, x_orbit_store(2,end)/1000, x_orbit_store(3,end)/1000, ...
      'rx', 'MarkerSize', 8, 'LineWidth', 2); 
xlabel('X_{ECI} [km]');
ylabel('Y_{ECI} [km]');
zlabel('Z_{ECI} [km]');
title('3D Satellite Orbital Path (ECI Frame)');
legend('Satellite Trajectory', 'Orbit Start', 'Orbit End', 'Location', 'best');
axis equal; 
view(3);  
hold off;

figure(11); % Position and Velocity
subplot(2,1,1); % pos
alt_store = vecnorm(x_orbit_store(1:3,:)) - 6378137;
plot(tspan, alt_store/1000, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Altitude [km]');
title('Orbital Altitude');

subplot(2,1,2); % vel
vel_store = vecnorm(x_orbit_store(4:6,:));
plot(tspan, vel_store, 'r-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Orbital Speed');

%% HCMG mode 
idx100 = tspan <= 100;
t100   = tspan(idx100);

figure(12); % HCMG Gimbal & Flywheel Actuation (zoomed)
tiledlayout(2,2,'TileSpacing','compact');

nexttile; % Gimbal rate
plot(t100, rad2deg(delta_dot_store_actual(1,idx100)), 'r-', ...
     t100, rad2deg(delta_dot_store_actual(2,idx100)), 'g-', ...
     t100, rad2deg(delta_dot_store_actual(3,idx100)), 'b-', ...
     t100, rad2deg(delta_dot_store_actual(4,idx100)), 'k-', 'LineWidth', lw);
grid on; ylabel('Gimbal Rate [deg/s]');
legend('\delta_1','\delta_2','\delta_3','\delta_4');

nexttile; % Gimbal angle
plot(t100, rad2deg(delta_store(1,idx100)), 'r-', ...
     t100, rad2deg(delta_store(2,idx100)), 'g-', ...
     t100, rad2deg(delta_store(3,idx100)), 'b-', ...
     t100, rad2deg(delta_store(4,idx100)), 'k-', 'LineWidth', lw);
grid on; ylabel('Gimbal Angle [deg]');
legend('\delta_1','\delta_2','\delta_3','\delta_4');

nexttile; % Flywheel accel
plot(t100, Omega_dot_store_actual(1,idx100), 'r-', ...
     t100, Omega_dot_store_actual(2,idx100), 'g-', ...
     t100, Omega_dot_store_actual(3,idx100), 'b-', ...
     t100, Omega_dot_store_actual(4,idx100), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Flywheel Acceleration [rad/s^2]');
legend('\Omega_1','\Omega_2','\Omega_3','\Omega_4');

nexttile; % Flywheel speed
plot(t100, Omega_store(1,idx100), 'r-', ...
     t100, Omega_store(2,idx100), 'g-', ...
     t100, Omega_store(3,idx100), 'b-', ...
     t100, Omega_store(4,idx100), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('Flywheel Speed [rad/s]');
legend('\Omega_1','\Omega_2','\Omega_3','\Omega_4');

sgtitle('HCMG Gimbal & Flywheel Actuation (First 100 s)');

figure(13); % HCMG Attitude Error, Torque, Singularity & Mode (zoomed)
tiledlayout(4,1,'TileSpacing','compact');

nexttile; % Attitude error quaternion
plot(t100, Ee_store(1,idx100), 'r-', t100, Ee_store(2,idx100), 'g-', ...
     t100, Ee_store(3,idx100), 'b-', t100, ne_store(idx100), 'k-', 'LineWidth', lw);
grid on; ylabel('Error Quaternion Component [-]');
legend('\epsilon_{e1}','\epsilon_{e2}','\epsilon_{e3}','\eta_e');

nexttile; % Torque
plot(t100, tau_cmd_store(1,idx100)*1000, 'r--', ...
     t100, tau_cmd_store(2,idx100)*1000, 'g--', ...
     t100, tau_cmd_store(3,idx100)*1000, 'b--', ...
     t100, tau_store(1,idx100)*1000,     'r-',  ...
     t100, tau_store(2,idx100)*1000,     'g-',  ...
     t100, tau_store(3,idx100)*1000,     'b-',  'LineWidth', lw);
grid on; ylabel('Torque [mNm]');
title('Commanded vs Actual HCMG Torque');
legend('\tau_{cmd,x}','\tau_{cmd,y}','\tau_{cmd,z}', ...
       '\tau_{act,x}','\tau_{act,y}','\tau_{act,z}');

nexttile; % Singularity parameters
plot(t100, S_CMG_store(idx100), 'b-', t100, S_RW_store(idx100), 'r-', 'LineWidth', lw);
grid on; ylabel('Singularity Parameter [-]');
title('Singularity Parameters');
legend('S_{CMG}','S_{RW}');

nexttile; % Mode switch
plot(t100, alpha_store(idx100), 'k-', 'LineWidth', lw);
grid on; xlabel('Time [s]'); ylabel('\alpha [-]');
title('Mode Switch Parameter \alpha');
ylim([-0.1 1.1]);

sgtitle('HCMG Attitude Error, Torque & Singularity (First 100 s)');