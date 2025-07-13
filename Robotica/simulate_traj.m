function [tgrid,q,dq,ddq,tau,e,q_ds,dq_ds,ddq_ds] = simulate_traj( ...
            robot,M_sym,C_sym,G_sym, Kp,Kd,Ki, coeffs, tf, N)
%SIMULATE_TRAJ  Integra a dinâmica do manipulador com PID clássico
%               usando método de Euler explícito.

    n     = robot.dof;
    tgrid = linspace(0,tf,N);
    dt    = tgrid(2)-tgrid(1);

    %% ---------- trajetória desejada (function handles) ----------
    syms t real
    [q_sym,dq_sym,ddq_sym] = calc_traj(coeffs,t);
    qd_fun   = matlabFunction(q_sym,   'Vars', t);
    dqd_fun  = matlabFunction(dq_sym,  'Vars', t);
    ddqd_fun = matlabFunction(ddq_sym, 'Vars', t);

    q_ds   = qd_fun(tgrid).';
    dq_ds  = dqd_fun(tgrid).';
    ddq_ds = ddqd_fun(tgrid).';

    %% ---------- cria M,C,G(q,dq) numericamente ------------------
    qvec  = sym('q',[n 1],'real');
    dqvec = sym('dq',[n 1],'real');

    M_fun = matlabFunction(subs(M_sym,  robot.q,             qvec),  'Vars',{qvec});
    C_fun = matlabFunction(subs(C_sym, [robot.q;robot.dq], [qvec;dqvec]), 'Vars',{qvec,dqvec});
    G_fun = matlabFunction(subs(G_sym,  robot.q,             qvec),  'Vars',{qvec});

    %% ---------- buffers -----------------------------------------
    q   = zeros(N,n); q(1,:) = q_ds(1,:);
    dq  = zeros(N,n);
    ddq = zeros(N,n);
    tau = zeros(N,n);
    e   = zeros(N,n);
    e_int = zeros(1,n);

    %% ---------- laço Euler --------------------------------------
    for k = 1:N-1
        %% erro e PID -------------------------------------------
        e(k,:)  = q_ds(k,:)  - q(k,:);
        edot    = dq_ds(k,:) - dq(k,:);
        e_int   = e_int + e(k,:)*dt;

        %% dinâmica no passo k ----------------------------------
        M = M_fun(q(k,:).');               % n×n
        C = C_fun(q(k,:).',dq(k,:).');     % n×n
        G = G_fun(q(k,:).').';             % 1×n row

        tau(k,:)   = (Kp*e(k,:).' + Kd*edot.' + Ki*e_int.' + G.').';
        rhs        = tau(k,:).' - C*dq(k,:).';

        % inversão robusta
        if rcond(M) < 1e-10
            ddq(k+1,:) = (pinv(M)*rhs).';
        else
            ddq(k+1,:) = (M\rhs).';
        end

        %% integração Euler ------------------------------------
        dq(k+1,:) = dq(k,:) + ddq(k+1,:)*dt;
        q(k+1,:)  = q(k,:)  + dq(k+1,:)*dt;
    end

    e(end,:) = q_ds(end,:) - q(end,:);
end

