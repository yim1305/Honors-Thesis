function f = mean2true(M, e)
    % Solve Kepler's equation for E
    % Use Newton-Raphson
    E = M;  % initial guess
    for k = 1:100
        E_new = E - (E - e*sin(E) - M) / (1 - e*cos(E));
        if abs(E_new - E) < 1e-12
            break;
        end
        E = E_new;
    end
    
    % Eccentric anomaly to true anomaly
    f = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
end