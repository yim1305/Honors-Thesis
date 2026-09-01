function euler = quat2eul(q)    
    Ex = q(1,:);
    Ey = q(2,:);
    Ez = q(3,:);
    n  = q(4,:);
    roll  = atan2(2*(n.*Ex + Ey.*Ez), 1 - 2*(Ex.^2 + Ey.^2));
    pitch = asin(max(-1, min(1, 2*(n.*Ey - Ez.*Ex))));
    yaw   = atan2(2*(n.*Ez + Ex.*Ey), 1 - 2*(Ey.^2 + Ez.^2));
    euler = [roll; pitch; yaw]; 
end