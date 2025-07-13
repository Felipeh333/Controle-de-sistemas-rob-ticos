function [q,qd,qdd] = calc_traj(coeffs,t)
%CALC_TRAJ  Avalia posição, velocidade e aceleração no instante t.
%   coeffs : 6×n, saída de coeff_traj
%   t      : escalar, vetor ou sym
%   q,qd,qdd : n×size(t) (mesma forma de t)

    t = sym(t);                   % permite usar tanto numérico quanto sym
    a0 = coeffs(1,:).'; a1 = coeffs(2,:).'; a2 = coeffs(3,:).';
    a3 = coeffs(4,:).'; a4 = coeffs(5,:).'; a5 = coeffs(6,:).';

    q   = a0 + a1.*t + a2.*t.^2 + a3.*t.^3 + a4.*t.^4 + a5.*t.^5;
    qd  = a1 + 2*a2.*t + 3*a3.*t.^2 + 4*a4.*t.^3 + 5*a5.*t.^4;
    qdd = 2*a2 + 6*a3.*t + 12*a4.*t.^2 + 20*a5.*t.^3;
end