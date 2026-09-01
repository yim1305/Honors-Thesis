function [qd, wd, wd_dot, Rd] = compute_nadir_reference(r_ECI, v_ECI)

    %% Position/velocity magnitudes

    r_mag = norm(r_ECI);

    %% LVLH frame axes

    r_hat = r_ECI / r_mag;

    h_vec = cross(r_ECI, v_ECI);
    h_mag = norm(h_vec);
    h_hat = h_vec / h_mag;

    % Desired body axes expressed in ECI

    zd = -r_hat;          % nadir
    yd = -h_hat;          % opposite orbit normal
    xd = cross(yd, zd);   % along-track

    xd = xd / norm(xd);

    %% Desired rotation matrix

    % Columns = body axes expressed in ECI
    Rd = [xd yd zd];

    %% Rotation matrix -> quaternion
    % Quaternion format: [Ex; Ey; Ez; n]

    qd = rotm2quat_scalar_last(Rd);

    %% Desired angular velocity (LVLH)

    omega_mag = h_mag / r_mag^2;

    wd = [0;
         -omega_mag;
          0];

    %% Desired angular acceleration (LVLH)

    r_dot = dot(r_ECI, v_ECI) / r_mag;

    omega_dot_mag = -2*(r_dot/r_mag)*omega_mag;

    wd_dot = [0;
             -omega_dot_mag;
              0];

end