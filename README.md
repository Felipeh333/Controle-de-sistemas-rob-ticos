# MATLAB Robot Dynamics & Control Toolkit

Scripts & funções principais:

- `Robot.m`           : classe – cinemática/dinâmica simbólica.
- `trajectory.m`      : `coeff_traj` + `calc_traj` (polinômio 5ª ordem).
- `dyn_on_traj.m`     : avalia `M`, `C`, `G` ao longo de `tgrid`.
- `control_gain.m`    : calcula `Kp`, `Kd`, `Ki` a partir de ζ e Tₛ.
- `simulate_traj.m`   : integra dinâmica com PID clássico (Euler).
- `plo_results.m`     : gera figuras (q, dq, ddq, erro, τ).
- `main.m`            : pipeline completo.

## Requisitos
    MATLAB R2023a + Symbolic Math Toolbox

## Execução rápida
    >> main   % roda tudo e plota resultados

## Fluxo resumido
    [Modelo simbólico] -> [Parâmetros numéricos] -> [Trajetória 5ª] ->
    [M,C,G(t)] -> [Ganhos PID] -> [Simulação] -> [Plots]
