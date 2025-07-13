function coeffs = coeff_traj(q0,qf,tf)
%COEFF_TRAJ  Coeficientes do polinômio de 5ª ordem para cada junta.
%   q0, qf : vetores coluna (n×1)
%   tf     : tempo total [s]
%   coeffs : matriz 6×n com [a0;…;a5] por coluna

    D = qf - q0;
    n = length(q0);
    coeffs = zeros(6,n);

    coeffs(1,:) = q0.';           % a0
    coeffs(2,:) = 0;              % a1
    coeffs(3,:) = 0;              % a2
    coeffs(4,:) = 10*D.'/tf^3;    % a3
    coeffs(5,:) = -15*D.'/tf^4;   % a4
    coeffs(6,:) = 6*D.'/tf^5;     % a5
end