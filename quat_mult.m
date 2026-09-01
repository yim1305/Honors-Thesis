function q = quat_mult(q1, q2)
    e1 = q1(1:3);
    n1 = q1(4);
    e2 = q2(1:3);
    n2 = q2(4);
    e = n1*e2 + n2*e1 + cross(e1,e2);
    n = n1*n2 - dot(e1,e2);
    q = [e; n];
    q = q / norm(q);
end
