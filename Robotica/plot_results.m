function plo_results(t,q,dq,ddq,qd,dqd,ddqd,e,tau)
%PLO_RESULTS  Cria figuras padrão (posições, velocidades, acelerações,
%erro e torques) para um manipulador com n juntas.

    n = size(q,2);
    jointNames = arrayfun(@(i) sprintf('q_{%d}',i), 1:n, 'UniformOutput',false);

    % unidades (replicadas se n>3)
    U_q   = repmat({'[rad]','[rad]','[m]'}, 1, ceil(n/3));  U_q   = U_q(1:n);
    U_dq  = repmat({'[rad/s]','[rad/s]','[m/s]'}, 1, ceil(n/3)); U_dq  = U_dq(1:n);
    U_ddq = repmat({'[rad/s^2]','[rad/s^2]','[m/s^2]'}, 1, ceil(n/3));U_ddq = U_ddq(1:n);

    plotGroup(q,  qd,  'Posições',       U_q  );
    plotGroup(dq, dqd, 'Velocidades',    U_dq );
    plotGroup(ddq,ddqd,'Acelerações',    U_ddq);
    plotGroup(e,  [],  'Erro de posição',U_q  );

    figure('Name','Torques','NumberTitle','off');
    plot(t, tau, 'LineWidth',1.2);
    grid on; xlabel('Tempo [s]'); ylabel('Torque / Força');
    title('Torques / força aplicada');
    legend(jointNames,'Location','best');

    %% ----------  função interna ----------------------------------
    function plotGroup(realData,refData,titleStr,units)
        figure('Name',titleStr,'NumberTitle','off');
        for i = 1:n
            subplot(n,1,i); hold on;
            plot(t, realData(:,i), 'LineWidth',1.4);
            if ~isempty(refData)
                plot(t, refData(:,i), '--k','LineWidth',0.8);
                legend({'real','ref'},'Location','best');
            end
            grid on;
            xlabel('Tempo [s]');
            ylabel(sprintf('%s %s', jointNames{i}, units{i}));
            if i==1, title(titleStr); end
        end
    end
end


