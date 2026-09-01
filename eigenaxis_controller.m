function [torque, Ee, ne, Re_bd] = eigenaxis_controller(q, qd, w, wd, wd_dot, J, Kp, Kd, h)
    E  = q(1:3);    n  = q(4);
    Ed = qd(1:3);   nd = qd(4);
    Matrix = [nd*eye(3) - skew(Ed), -Ed;
              Ed',                   nd];
    result = Matrix * [E; n];     
    Ee = result(1:3);
    ne = result(4);
    if ne < 0                    
        Ee = -Ee;
        ne = -ne;
    end
    Re_bd  = quat2rotm_body(q)' * quat2rotm_body(qd);  
    we     = w - Re_bd * wd;
    wd_dot = Re_bd * wd_dot;
 
    torque = -skew(w)*(J*w + h) - J*wd_dot + Kp*J*Ee + Kd*J*we;
end