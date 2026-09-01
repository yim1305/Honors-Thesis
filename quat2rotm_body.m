function R = quat2rotm_body(q)

% body frame to ECI

    Ex = q(1);
    Ey = q(2);
    Ez = q(3);
    n = q(4);

    R = [(1-2*(Ey^2+Ez^2)),  2*(Ex*Ey-n*Ez),   2*(Ex*Ez+n*Ey);
          2*(Ex*Ey+n*Ez),   (1-2*(Ex^2+Ez^2)),  2*(Ey*Ez-n*Ex);
          2*(Ex*Ez-n*Ey),    2*(Ey*Ez+n*Ex),   (1-2*(Ex^2+Ey^2))];
end