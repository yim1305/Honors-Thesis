function q = eul2quat(roll, pitch, yaw)

% roll - x, pitch - y, yaw - z

    % Compute half-angles
    cr = cos(roll / 2);
    sr = sin(roll / 2);
    cp = cos(pitch / 2);
    sp = sin(pitch / 2);
    cy = cos(yaw / 2);
    sy = sin(yaw / 2);

    % Calculate quaternion components ZYX sequence
    qx = sr * cp * cy - cr * sp * sy;
    qy = cr * sp * cy + sr * cp * sy;
    qz = cr * cp * sy - sr * sp * cy;
    qw = cr * cp * cy + sr * sp * sy;

    q = [qx; qy; qz; qw];
    
    % normalize
    q = q / norm(q);
end