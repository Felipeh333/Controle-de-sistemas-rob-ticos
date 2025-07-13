function [Kp,Kd,Ki,wn,p] = control_gain(pid,zeta,ts,M_max,C_max)
%CONTROL_GAIN  Calcula matrizes de ganho PID (clássico ou torque calculado)
%   pid  : 'classico' ou 'tc'
%   zeta : razão de amortecimento desejada
%   ts   : settling time (5%)
%   M_max, C_max : matrizes dinâmicas nos piores casos

    wn = 3/(zeta*ts);
    p  = 4*wn;

    switch lower(pid)
        case 'classico'
            Kp = diag(diag((wn^2 + 2*zeta*wn*p) * M_max));
            Kd = diag(diag((2*zeta*wn + p)   * M_max - C_max));
            Ki = diag(diag((wn^2*p) * M_max));
        case 'tc'
            n  = size(M_max,1);
            Kp = diag( repmat(wn^2 + 2*zeta*wn*p, 1, n) );
            Kd = diag( repmat(2*zeta*wn + p,     1, n) );
            Ki = diag( repmat(wn^2*p,             1, n) );
        otherwise
            error('pid deve ser "classico" ou "tc"');
    end
end