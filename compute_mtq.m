function [tau_mtq, m_dipole] = compute_mtq(h_excess, sat, B_body)
    B_mag2   = dot(B_body, B_body);
    if B_mag2 < 1e-20
        tau_mtq  = zeros(3,1);
        m_dipole = zeros(3,1);
        return;
    end
    m_dipole = sat.mtq.k_desat * cross(h_excess, B_body) / B_mag2;
    % clamp
    for k = 1:3
        m_dipole(k) = max(-sat.mtq.m_max, min(sat.mtq.m_max, m_dipole(k))); % Am²
    end
    % torque
    tau_mtq = cross(m_dipole, B_body); % Nm
end