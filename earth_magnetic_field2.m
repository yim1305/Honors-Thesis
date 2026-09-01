function B_ECI = earth_magnetic_field2(r_ECI, t)
    mu_mag      = 7.94e22;       
    mu0         = 4*pi*1e-7;      
    gamma       = deg2rad(11.3);
    omega_earth = 7.2921159e-5;   
    m_hat_ECEF = - [sin(gamma); 0; cos(gamma)];
    theta = omega_earth * t;    
    R_z = [cos(theta) -sin(theta) 0;
           sin(theta)  cos(theta) 0;
           0           0          1];
    m_hat = R_z * m_hat_ECEF;
    r_mag = norm(r_ECI);
    r_hat = r_ECI / r_mag;
    B_ECI = (mu0*mu_mag/(4*pi*r_mag^3)) * ...
            (3*dot(m_hat, r_hat)*r_hat - m_hat);
end