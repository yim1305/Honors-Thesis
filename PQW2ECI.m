function x_ECI = PQW2ECI(x_PQW, parameters)
    i     = parameters(3);
    Omega = parameters(4);
    w     = parameters(5);
    C = [(cos(Omega)*cos(w) - sin(Omega)*sin(w)*cos(i)), (-cos(Omega)*sin(w) - sin(Omega)*cos(w)*cos(i)),  (sin(Omega)*sin(i));
         (sin(Omega)*cos(w) + cos(Omega)*sin(w)*cos(i)), (-sin(Omega)*sin(w) + cos(Omega)*cos(w)*cos(i)), (-cos(Omega)*sin(i));
         (sin(w)*sin(i)),                                (cos(w)*sin(i)),                                  (cos(i))];
    x_ECI = C * x_PQW;
end
