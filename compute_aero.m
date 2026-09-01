function [F_drag_body, a_drag_ECI, tau_drag_body] = compute_aero(q_curr, v_ECI, r_ECI, sat)

    %% Altitude and density
    R_earth = 6378137;
    h_alt   = norm(r_ECI) - R_earth;
    rho     = atmosphere_density(h_alt);

    %% Velocity in body frame
    R       = quat2rotm_body(q_curr);   % ECI to body 
    v_body  = R' * v_ECI;               % body
    v_mag   = norm(v_body);

    if v_mag < 1e-6
        F_drag_body = zeros(3,1);
        tau_drag_body    = zeros(3,1);
        return;
    end
    v_hat = v_body / v_mag;  

    %% Sum over all faces
    F_drag_body = zeros(3,1);
    tau_drag_body    = zeros(3,1);

    for k = 1:6
        n_hat     = sat.faces(k).normal;
        A_face    = sat.faces(k).area;
        cop       = sat.faces(k).cop - sat.CoM;  

        % Only exposed faces 
        cos_theta = dot(v_hat, n_hat); 
        if cos_theta <= 0
            continue;
        end

        % Drag force 
        dF = -0.5 * rho * sat.Cd * A_face * v_mag^2 * cos_theta * v_hat;

        F_drag_body = F_drag_body + dF;
        tau_drag_body    = tau_drag_body + cross(cop, dF);
    end

    F_drag_ECI = R * F_drag_body;    % body to ECI
    a_drag_ECI = F_drag_ECI / sat.m; 
end