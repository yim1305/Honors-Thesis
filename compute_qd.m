function [qd, wd, wd_dot, Rd] = compute_qd(r_ECI, v_ECI)

    %% Position/velocity magnitudes

    r_mag = norm(r_ECI);

    %% LVLH frame axes

    r_hat = r_ECI / r_mag;

    h_vec = cross(r_ECI, v_ECI);
    h_mag = norm(h_vec);
    h_hat = h_vec / h_mag;

    %% Desired body axes expressed in ECI
    yd = -r_hat;          % nadir (now on Y)
    zd = -h_hat;           % opposite orbit normal (now on Z)
    xd = cross(yd, zd);    % completes right-handed frame
    xd = xd / norm(xd);
    
    %% Desired rotation matrix
    Rd = [xd yd zd];       % columns unchanged: [X Y Z] body axes in ECI


    %% Rotation matrix -> quaternion
    % Quaternion format: [Ex; Ey; Ez; n]

    qd = rotm2quat_scalar_last(Rd);

    %% Desired angular velocity (LVLH)

    omega_mag = h_mag / r_mag^2;

    %% Desired angular velocity (now in Z slot)
    wd = [0;
          0;
         -omega_mag];

    %% Desired angular acceleration (LVLH)

    r_dot = dot(r_ECI, v_ECI) / r_mag;

    omega_dot_mag = -2*(r_dot/r_mag)*omega_mag;

    %% Desired angular acceleration
    wd_dot = [0;
              0;
             -omega_dot_mag];

end