function r_I_B = ECI2BODY(q, nadir_I)
    qx = q(1);
    qy = q(2);
    qz = q(3);
    qw = q(4);
    R = [1 - 2*(qy^2 + qz^2), 2*(qx*qy + qw*qz),   2*(qx*qz - qw*qy);
         2*(qx*qy - qw*qz),   1 - 2*(qx^2 + qz^2), 2*(qy*qz + qw*qx);
         2*(qx*qz + qw*qy),   2*(qy*qz - qw*qx),   1 - 2*(qx^2 + qy^2)];
    r_I_B = R * nadir_I;
end
