function [M_series,C_series,G_series,M_max,C_max] = dyn_on_traj(robot,M_sym,C_sym,G_sym,coeffs,tgrid)
% Avalia M(q), C(q,qdot) e G(q) ao longo de uma trajetória polinomial
% de 5ª ordem e devolve as séries + o maior valor (norma Frobenius).

    syms t real

    % ----- trajetória simbólica -----
    [q_d, dq_d, ddq_d] = calc_traj(coeffs, t);

    % ----- pares {sym, val} para evalSubs -----
    varsSub = [robot.q ; robot.dq ; robot.ddq];
    valsSub = [q_d     ; dq_d     ; ddq_d   ];
    pairsQC = [num2cell(varsSub), num2cell(valsSub)];

    % ----- expressões só em função de t -----
    M_t = robot.evalSubs(M_sym, pairsQC);
    C_t = robot.evalSubs(C_sym, pairsQC);
    G_t = robot.evalSubs(G_sym, pairsQC);

    % ----- function handles numéricas -----
    M_fun = matlabFunction(M_t, 'Vars', t);
    C_fun = matlabFunction(C_t, 'Vars', t);
    G_fun = matlabFunction(G_t, 'Vars', t);

    N = numel(tgrid);
    n = robot.dof;
    M_series = zeros(n,n,N);
    C_series = zeros(n,n,N);
    G_series = zeros(n,1,N);

    for k = 1:N
        tk = tgrid(k);
        M_series(:,:,k) = M_fun(tk);
        C_series(:,:,k) = C_fun(tk);
        G_series(:,:,k) = G_fun(tk);
    end

    % ----- normas e máximos -----
    M_norm = arrayfun(@(idx) norm(M_series(:,:,idx),'fro'), 1:N);
    C_norm = arrayfun(@(idx) norm(C_series(:,:,idx),'fro'), 1:N);

    [~,idxM] = max(M_norm);
    [~,idxC] = max(C_norm);

    M_max = M_series(:,:,idxM);
    C_max = C_series(:,:,idxC);
end

