function A = matA(theta, delta)
    delta1 = delta(1);
    delta2 = delta(2);
    delta3 = delta(3);
    delta4 = delta(4);
    A = [-cos(theta)*sin(delta1), -cos(delta2),             cos(theta)*sin(delta3),  cos(delta4);
          cos(delta1),           -cos(theta)*sin(delta2),  -cos(delta3),             cos(theta)*sin(delta4);
          sin(theta)*sin(delta1), sin(theta)*sin(delta2),   sin(theta)*sin(delta3),  sin(theta)*sin(delta4)];
end


%% Original
%{
function A = matA(theta, delta)

    delta1 = delta(1);
    delta2 = delta(2);
    delta3 = delta(3);
    delta4 = delta(4);

    A = [-cos(theta)*sin(delta1), cos(delta2),             cos(theta)*sin(delta3), -cos(delta4);
          cos(delta1),           -cos(theta)*sin(delta2), -cos(delta3),             cos(theta)*sin(delta4);
          sin(theta)*sin(delta1), sin(theta)*sin(delta2),  sin(theta)*sin(delta3),  sin(theta)*sin(delta4)];

end
%}