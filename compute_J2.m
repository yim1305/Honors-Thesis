function a_J2 = compute_J2(x)
    mu      = 3.986004418e14;
    J2      = 1.08263e-3;
    R_earth = 6378137;
    r       = x(1:3);
    z       = x(3);
    r_mag   = norm(r);
    factor  = -1.5 * J2 * mu * R_earth^2 / r_mag^4;
    a_J2 = factor * [(1 - 5*(z/r_mag)^2) * r(1)/r_mag;
                     (1 - 5*(z/r_mag)^2) * r(2)/r_mag;
                     (3 - 5*(z/r_mag)^2) * r(3)/r_mag];
end