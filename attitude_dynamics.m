function x_dot = attitude_dynamics(x, I, tau, h)
    E = [x(1); x(2); x(3)];
    n = x(4);
    w = [x(5); x(6); x(7)];
    E_dot = (0.5 * skew(E) * w) + (0.5 * n * w);
    n_dot = (-0.5 * w' * E);
    w_dot = I \ (-tau - skew(w)*(I*w + h));
    x_dot = [E_dot; n_dot; w_dot];
end