classdef Robot
    % ROBOT – Dinâmica simbólica de um manipulador (ex.: PRR/RRP)
    %   Portado de Python/SymPy para MATLAB Symbolic Math Toolbox.
    
    properties
        % ------------ configurações gerais -------------
        name          string = "Cookboat"
        jointTypes    char        % ex.: 'PRR'
        dof           double      % número de juntas
        
        % ------------ parâmetros DH --------------------
        thetas   sym
        ds       sym
        a        sym
        alphas   sym
        
        % ------------ variáveis de estado --------------
        t    sym
        q    sym
        dq   sym
        ddq  sym
        
        % ------------ parâmetros dinâmicos -------------
        masses        sym          % [dof×1]
        r_cis_local   cell         % {dof} cada 3×1
        inertias      cell         % {dof} cada 3×3
        gVec          sym          % 3×1
        
        % ------------ resultados -----------------------
        A_i   cell   % DH individuais
        T_0i  cell   % base→i
        J     sym    % Jacobiano EE
        J_ci  cell   % Jacobianos por COM
        M     sym    % matriz de inércia
        C     sym    % matriz de Coriólis
        G     sym    % vetor de gravidade
        tau   sym    % τ = M*q̈ + C*q̇ + G
    end
    
    methods
        %==============================================================
        function obj = Robot(configuration, dh_table, varargin)
            % Construtor:
            %   Robot("PRR", dh_table, 'masses', m, 'r_cis', {...}, …)
            
            % --- configuração básica ---
            obj.jointTypes = char(configuration);
            obj.dof        = length(obj.jointTypes);
            obj.t          = sym('t','real');
            
            obj.q   = sym('q',  [obj.dof 1], 'real');
            obj.dq  = sym('dq', [obj.dof 1], 'real');
            obj.ddq = sym('ddq',[obj.dof 1], 'real');
            
            % --- valores-padrão dinâmicos ---
            obj.masses      = sym('m', [obj.dof 1], 'real');
            obj.r_cis_local = arrayfun(@(~) sym(zeros(3,1)), 1:obj.dof, ...
                                       'UniformOutput', false);
            obj.inertias    = arrayfun(@(~) sym('I', [3 3], 'real'), ...
                                       1:obj.dof, 'UniformOutput', false);
            obj.gVec        = [0;0;-sym('g','real')];
            
            % --- parse pares Nome–Valor ---
            if mod(numel(varargin),2) ~= 0
                error("Parâmetros devem vir em pares nome-valor.");
            end
            for k = 1:2:numel(varargin)
                name  = lower(string(varargin{k}));
                value = varargin{k+1};
                switch name
                    case "masses",       obj.masses      = value;
                    case {"r_cis","r_cis_local"}
                                         obj.r_cis_local = value;
                    case "inertias",     obj.inertias    = value;
                    case "g_vec",        obj.gVec        = value;
                    otherwise
                        error("Parâmetro desconhecido: %s", name);
                end
            end
            
            % --- parâmetros DH ---
            obj = obj.assignDH(dh_table);
            
            % --- cinemática & dinâmica ---
            obj = obj.forwardKinematics();
            obj = obj.computeJacobians();
            obj = obj.computeDynamics();
        end
        %--------------------------------------------------------------
        function obj = assignDH(obj, dh)
            obj.thetas = dh(:,1);
            obj.ds     = dh(:,2);
            obj.a      = dh(:,3);
            obj.alphas = dh(:,4);
            
            for i = 1:obj.dof
                qi = symfun( sym(sprintf('q%d',i),'real'), obj.t );
                if obj.jointTypes(i)=='R'
                    obj.thetas(i) = qi;
                else % 'P'
                    obj.ds(i) = qi;
                end
            end
        end
        %--------------------------------------------------------------
        function A = dhMatrix(~, th, d, a, al)
            A = [ cos(th) -sin(th)*cos(al)  sin(th)*sin(al)  a*cos(th)
                  sin(th)  cos(th)*cos(al) -cos(th)*sin(al)  a*sin(th)
                  0        sin(al)          cos(al)          d
                  0        0                0                1 ];
        end
        %--------------------------------------------------------------
        function obj = forwardKinematics(obj)
            obj.A_i  = cell(obj.dof,1);
            obj.T_0i = cell(obj.dof,1);
            T = sym(eye(4));
            for i = 1:obj.dof
                obj.A_i{i} = obj.dhMatrix(obj.thetas(i),obj.ds(i), ...
                                          obj.a(i),obj.alphas(i));
                T          = T * obj.A_i{i};
                obj.T_0i{i}= T;
            end
        end
        %--------------------------------------------------------------
        function z = zi(obj,i)
            if i==0,  z = [0;0;1];
            else,     z = obj.T_0i{i}(1:3,3);
            end
        end
        function p = pi(obj,i)
            if i==0,  p = sym(zeros(3,1));
            else,     p = obj.T_0i{i}(1:3,4);
            end
        end
        %--------------------------------------------------------------
        function obj = computeJacobians(obj)
            pe = obj.T_0i{end}(1:3,4);
            J  = sym(zeros(6,obj.dof));
            for i = 1:obj.dof
                if obj.jointTypes(i)=='R'
                    Jv = cross(obj.zi(i-1), pe - obj.pi(i-1));
                    Jw = obj.zi(i-1);
                else
                    Jv = obj.zi(i-1);
                    Jw = sym(zeros(3,1));
                end
                J(:,i) = [Jv;Jw];
            end
            obj.J = simplify(J);
            
            % Jacobianos em cada centro de massa
            obj.J_ci = cell(obj.dof,1);
            for link = 1:obj.dof
                pc = obj.T_0i{link} * [obj.r_cis_local{link};1];
                pc = pc(1:3);
                Jc = sym(zeros(6,obj.dof));
                for i = 1:obj.dof
                    if i>link && obj.jointTypes(i)=='P'
                        continue
                    end
                    if obj.jointTypes(i)=='R'
                        Jv = cross(obj.zi(i-1), pc - obj.pi(i-1));
                        Jw = obj.zi(i-1);
                    else
                        Jv = obj.zi(i-1);
                        Jw = sym(zeros(3,1));
                    end
                    Jc(:,i) = [Jv;Jw];
                end
                obj.J_ci{link} = simplify(Jc);
            end
        end
        %--------------------------------------------------------------
        function obj = computeDynamics(obj)
            % ----- Inércia M(q) -----
            M = sym.zeros(obj.dof);
            for i = 1:obj.dof
                Jv = obj.J_ci{i}(1:3,:);
                Jw = obj.J_ci{i}(4:6,:);
                M  = M + obj.masses(i)*(Jv.'*Jv) + ...
                          (Jw.' * obj.inertias{i} * Jw);
            end
            obj.M = simplify(M);
            
            % ----- Coriólis C(q,q̇) -----
            n = obj.dof;
            C = sym.zeros(n);
            for i = 1:n
                for j = 1:n
                    cij = sym(0);
                    for k = 1:n
                        c = 0.5*( diff(M(i,k),obj.q(j)) + ...
                                  diff(M(i,j),obj.q(k)) - ...
                                  diff(M(j,k),obj.q(i)) );
                        cij = cij + c*obj.dq(k);
                    end
                    C(i,j) = cij;
                end
            end
            obj.C = simplify(C);
            
            % ----- Gravidade G(q) -----
            U = sym(0);
            for i = 1:obj.dof
                pc = obj.T_0i{i} * [obj.r_cis_local{i};1];
                U  = U + obj.masses(i) * dot(obj.gVec, pc(1:3));
            end
            obj.G = simplify( jacobian(U, obj.q).' );
            
            % ----- τ = M*q̈ + C*q̇ + G -----
            obj.tau = simplify( obj.M*obj.ddq + obj.C*obj.dq + obj.G );
        end
    end
    
    %===================== utilitários ===============================
    methods
        function exprNum = evalSubs(~, expr, pairs)
            % Substitui variáveis simbólicas por valores numéricos.
            %
            %   • pairs pode ser:
            %         – struct:   s.m1 = 0.9;  s.g = 9.81; …
            %         – célula:   {sym  val;   sym  val;  …}
            
            % ---- struct → célula ------------------------------------
            if isstruct(pairs)
                names = fieldnames(pairs);
                pairs = cellfun(@(n){sym(n), pairs.(n)}, names, ...
                                'UniformOutput', false);
                pairs = vertcat(pairs{:});
            elseif ~iscell(pairs)
                error('pairs deve ser struct ou célula {sym val; …}');
            end
            
            % ---- separa em dois vetores (vars, vals) ---------------
            vars = [pairs{:,1}];              % 1×N sym
            vals = arrayfun(@sym, [pairs{:,2}]);  % 1×N sym
            
            % ---- substitui -----------------------------------------
            exprNum = simplify( subs(expr, vars, vals) );
        end
    end
end

