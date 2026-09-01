function B = matB(theta, delta)
    delta1 = delta(1);
    delta2 = delta(2);
    delta3 = delta(3);
    delta4 = delta(4);
    B = [ cos(theta)*cos(delta1), -sin(delta2),            -cos(theta)*cos(delta3),  sin(delta4);
          sin(delta1),             cos(theta)*cos(delta2), -sin(delta3),            -cos(theta)*cos(delta4);
         -sin(theta)*cos(delta1), -sin(theta)*cos(delta2), -sin(theta)*cos(delta3), -sin(theta)*cos(delta4)];
end
