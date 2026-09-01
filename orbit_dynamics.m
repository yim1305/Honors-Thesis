function xdot = orbit_dynamics(x, a_drag, a_J2)
    mu = 3.986004418e14;
    rx = x(1);
    ry = x(2);
    rz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);
    v = [vx; vy; vz];
    r = sqrt(rx^2 + ry^2 + rz^2);
    ax = -mu * rx / r^3;
    ay = -mu * ry / r^3;
    az = -mu * rz / r^3;
    a = [ax; ay; az];
    xdot = [v; a + a_drag + a_J2];
end
