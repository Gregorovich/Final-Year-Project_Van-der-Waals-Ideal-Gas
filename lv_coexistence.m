fprintf('Calculating\n')

% Water:
%  a = 5.536   [ L^2 bar mol^-1 ]
%  b = 0.03049 [ L mol^-1 ]

% Console Information
syms           Vr Pr Tr pr X Vr_g Vr_l
P_eqn_Vol      = Pr == -(8*Tr*Vr^2 - 9*Vr + 3) / (Vr^2 - 3*Vr^3);
P_eqn_Vol_3bic =  0 == 3*Pr*Vr^3 - (8*Tr+Pr)*Vr^2 + 9*Vr - 3 ;
% Pg = Pl; Gg = Gl
% P_dens       = 8*pr*Tr/(3   -pr) - 3*pr^2
% P_Vol        = 8*Tr   /(3*Vr-1 ) - 3/Vr^2
P_eqn_Vol_PG   = -(8*Tr*Vr_l^2 - 9*Vr_l + 3) / (Vr_l^2 - 3*Vr_l^3) ...
              == -(8*Tr*Vr_g^2 - 9*Vr_g + 3) / (Vr_g^2 - 3*Vr_g^3);
% G_dens       = -Tr*log[(3     -pr )/pr ] + Tr*p /(3     -pr ) - 9*pr /4        - Tr*log(X)
% G_Vol        = -Tr*log( 3*Vr  -1       ) + Tr   /(3*Vr  -1  ) - 9   /(4*Vr   ) - Tr*log(X)
G_eqn_Vol_PG   = -Tr*log( 3*Vr_g-1       ) + Tr   /(3*Vr_g-1  ) - 9   /(4*Vr_g ) - Tr*log(X) ...
              == -Tr*log( 3*Vr_l-1       ) + Tr   /(3*Vr_l-1  ) - 9   /(4*Vr_l ) - Tr*log(X);

fprintf('Solving this:\n%s\n%s\n\nAnd also solving Pg = Pl; Gg = Gl :\n            %s\n%s\n',P_eqn_Vol,P_eqn_Vol_3bic,P_eqn_Vol_PG,G_eqn_Vol_PG)

% Preallocation
t_end          = 10;
Vn             = zeros(3);
solns_Vr       = zeros(3,t_end^2);
solns_pr       = zeros(3,t_end^2);
solns_PG       = zeros(2,t_end);
Tr             = linspace(0,1,t_end);
update         = 1;

%% Solve Pressure-Volume-Temperature Cubic to find LV Coexistence Curve.
fprintf('\nP-V Cubic\n')
k = 1;
grph_end = 12; % multiple of 11
Vr_grph = linspace(0,1.1,grph_end);
pr_grph = linspace(0,1.1,grph_end);
Pr_grph = zeros(grph_end,1);
for i=0:t_end-1
    
    Tr = (i+1)/t_end;
    
    for j=1:t_end
        
        if rem(i,update)==0
            fprintf('.')
        end
        
        Pr = j/t_end;
        
        % 0 == 3*Vr^3 - Vr^2 - (8*Tr*Vr^2 - 9*Vr + 3)/Pr
        P_eqn_Vol = 0 == 3*Pr*Vr^3 -(8*Tr+Pr)*Vr^2 + 9*Vr -3 ;
        %P_eqn_Vol = [3*Pr*Vr^3 -(8*Tr+Pr)*Vr^2 9*Vr -3 ] ;
        
        Vn(1) = root(P_eqn_Vol,Vr,1);
        Vn(2) = root(P_eqn_Vol,Vr,2);
        Vn(3) = root(P_eqn_Vol,Vr,3);
        
        Vn = sort(Vn(:));
        
        V1 = Vn(1);
        V2 = Vn(2);
        V3 = Vn(3);
        
        solns_Vr(3*t_end*i+3*(j-1)+1) = V1;
        solns_Vr(3*t_end*i+3*(j-1)+2) = V2;
        solns_Vr(3*t_end*i+3*(j-1)+3) = V3;
        
        % 0 == 3*pr^3 - pr^2 - (8*Tr*pr^2 - 9*pr + 3)/Pr
        P_eqn_dens = 0 == 3*Pr*pr^3 -(8*Tr+Pr)*pr^2 + 9*pr -3 ;
        %P_eqn_dens = [3*Pr*pr^3 -(8*Tr+Pr)*pr^2 9*pr -3 ] ;
        
        pn(1) = root(P_eqn_dens,pr,1);
        pn(2) = root(P_eqn_dens,pr,2);
        pn(3) = root(P_eqn_dens,pr,3);
        
        pn = sort(pn(:));
        
        p1 = pn(1);
        p2 = pn(2);
        p3 = pn(3);
        
        solns_pr(3*t_end*i+3*(j-1)+1) = p1;
        solns_pr(3*t_end*i+3*(j-1)+2) = p2;
        solns_pr(3*t_end*i+3*(j-1)+3) = p3;
        
    end
end


% THIS BIT DDOESN'T QUITE WORK YET, WILL SOLVE ONCE THE GRAPH NO LONGER PLOTS STRAIGHT LINES
%% Solve Pg = Pl; Gg = Gl to find LV Coexistence Curve.
fprintf('\nAnd now Pg = Pl; Gg = Gl\n')
for k=1:1
    % Pg = Pl; Gg = Gl
    % P_dens        =  8 p  T/(3   -p ) - 3* p^2
    % P_Vol         =  8 T   /(3 V -1 ) - 3/ V^2
    % G_dens        =  -T  ln[(3   -p )/p ] + T p /(3   -p ) - 9p /4      - Tln(X)
    % G_Vol         =  -T  ln( 3 V -1     ) + T   /(3 V1-1 ) - 9  /(4 V ) - Tln(X)
    %   Tln(X) is a function of T produced by integrating P, and cancels out when
    %   solving the simultaneous equation, so is irrelavent.
    
    T = 1/t_end;
    options = optimset('Display','iter');
    
    F = { @(V) (8*T  /(3*V(1)-1) - 3/V(1)^2 - ( 8*T  /(3*V(2)-1) - 3/V(2)^2) ),
          @(V) ( -T*log( 3*V(1)-1   ) + T  /(3*V(2)-1) - 9 /(4*V(1)) - ( -T*log( 3*V(2)-1   ) + T  /(3*V(2)-1) - 9 /(4*V(2)) ) ) };
    
    -((4*N*Tr*a)/(9*b) - (N*a)/(3*Vr*b))/(3/Vr + (8*Tr*log(Vr - 1/3))/3 + (3*N*Vr*b*(8*Tr*Vr^2 - 9*Vr + 3))/(Vr^2 - 3*Vr^3))

    
    soln = fsolve(F, [ solns_Vr(1) , solns_Vr(3) ] , options);
    %soln = fsolve(F, [ solns_Vr(1) , solns_Vr(3) ] ) 
    
    
    
    solns_PG(1) = soln(1);
    solns_PG(2) = soln(2);

    for i=1:t_end-1
        
        if rem(i,update)==0
            fprintf('.')
        end
        
        T = (i+1)/t_end;
        
        F = { @(V) (8*T  /(3*V(1)-1) - 3/V(1)^2 - ( 8*T  /(3*V(2)-1) - 3/V(2)^2) ),
              @(V) ( -T*log( 3*V(1)-1   ) + T  /(3*V(1)-1) - 9 /(4*V(1)) - ( -T*log( 3*V(2)-1   ) + T  /(3*V(2)-1) - 9 /(4*V(2)) ) ) };
        
        soln = fsolve(F, [solns_PG(i*2),solns_PG(i*2+1)] , options);
        %soln = fsolve( F, [solns_PG(i*2),solns_PG(i*2+1)] )
        
        solns_PG(i*2+1) = soln(1);
        solns_PG(i*2+2) = soln(2);
        
    end
end


%% Plot Graphs
hold on

% Plot Pr/Vr Cubic:
Isotherms.m

% LV Coexistence from solving Pr/Vr cubic:
LV_P_cubic.m

% LV coexistence from Pg = Pl and Gg = Gl:
LV_PG_Smltns.m


hold off
fprintf('\n')
solns_trackd=solns_trackd'; % make it easier to view the data