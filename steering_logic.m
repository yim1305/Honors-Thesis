function [lambda, beta, delta_dot, Omega_dot] = steering_logic(t, lambda0, mu1, S_CMG, beta0, mu2, S_RW, Iw, Omega, alpha, A, B, tau, d)

    e1 = 0.001*sin(t+0.25);
    e2 = 0.001*sin(t+0.13);
    e3 = 0.001*sin(t+0.05);

    E = [1, e1, e2;
         e1, 1, e3;
         e2, e3, 1];

    lambda = lambda0*exp(-mu1 * S_CMG);
    beta   = beta0*exp(-mu2 * S_RW);

    Omega_diag = diag(Omega);

    cmg_torque_term = A' * ((A * A' + lambda * E) \ tau);
    term2           = cmg_torque_term + beta * d;
    delta_dot       = (alpha / Iw) * (Omega_diag \ term2);

    Omega_dot       = ((1 - alpha) / Iw) * (B' * ((B * B') \ tau));

end
