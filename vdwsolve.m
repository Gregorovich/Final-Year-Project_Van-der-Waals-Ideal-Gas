%function vdwsolve
% Rearranges the Van Der Waals Equation to remove a and b and make it
% independant of the gas. Also solves for Internal Energy, Enhalpy, and
% Gibbs Free Energy, and plots the graph of P vs ? and P vs V
%
% This is done by dP/dp = 0 and d2P/dp2 = 0 , and solving.
%
% Some fprintf's are used in place of comments here as I used excess 
% to act both as comments and to show the working in the
% command window when this script is run.


% Symbolic preallocation.
syms kB b a N T Tc Tr p pc pr V Vc Vr P Pc Pr E Ec Er Ec_eqn G Gc Gr A H Hc Hr

% Define Starting Equations:
P_dens  = P ==   p*kB*T * (1-b*p)^-1 -   a*p^2;
%dPdp   = 0 == - kB*T * (p-b)^-2 + 2*a*p^-3;
dPdp    = diff (P_dens,p);
%d2Pdp2 = 0 == 2*kB*T * (p-b)^-3 - 6*a*p^-4;
d2Pdp2  = diff (dPdp,p);

P_Vol   = P ==   N*kB*T * (V-N*b)^-1 -   N^2*a*V^-2;
%dPdV   = 0 == - N*kB*T * (V-N*b)^-2 + 2*N^2*a*V^-3;
dPdV    = diff (P_Vol,V);
%d2PdV2 = 0 == 2*N*kB*T * (V-N*b)^-3 - 6*N^2*a*V^-4;
d2PdV2  = diff (dPdV,V);
% MATLAB orders in preference of capitals and beginning of the alphabet
% first in subtractions, e.g.: (b-p) rather than (p-b), and (V-b), hence
% why these equations look different when MATLAB prints them. Haven't
% figured out why sometimes the first term is the V^-# term and sometimes
% the (V-b)^-# term (or equivalent density versions)

fprintf('Van der Waals Equation:\n          %s\t\t\t\t\t  %s\n  dP/dp = %s\t\t  dP/dV = %s\nd2P/dp2 = %s\td2P/dV2 = %s\n\n',P_dens,P_Vol,dPdp,dPdV,d2Pdp2,d2PdV2)


fprintf('Rearrange dPdp / dPdV and d2Pdp2 / d2PdV2 in terms of T\n')

dPdp_T   = solve (dPdp,  T);
d2Pdp2_T = solve (d2Pdp2,T);

dPdV_T   = solve (dPdV,  T);
d2PdV2_T = solve (d2PdV2,T);

fprintf('T = %s\tT = %s\nT = %s\tT = %s\n\n',dPdp_T,dPdV_T,d2Pdp2_T,d2PdV2_T)


% dPdp_T = d2Pdp2_T by definition - both are in terms of T
fprintf('Now find pc / Vc in terms of a and b\n')

smltns_eqn_dens = simplify(dPdp_T*(kB*b)/a)/(b*p-1)^2 == simplify(d2Pdp2_T*(kB*b)/a)/(b*p-1)^2;
smltns_eqn_Vol  = dPdV_T*(kB*V^4)*(V-N*b)^-2 == d2PdV2_T*(kB*V^4)*(V-N*b)^-2;
% Due to the presence of powers, MATLAB tries to solve this equation as
% either a quadratic or a cubic and as such finds multiple solutions,
% causing issues with the next formulae.
% Using "simplify (smltns_eqn)" does not work as the answer ends up being
% kB ~= 0 & p ~= 0 & (b == p | 3*b == p | a == 0)

pc_eqn = pc == solve (smltns_eqn_dens,p);
pc_val =      solve (smltns_eqn_dens,p);

Vc_eqn  = Vc == solve (smltns_eqn_Vol,V);
Vc_val  =       solve (smltns_eqn_Vol,V);

fprintf('%s\t%s\n   %s\t\t\t\t %s\n\n',smltns_eqn_dens,smltns_eqn_Vol,pc_eqn,Vc_eqn)


fprintf('Substitute pc into dPdp to find Tc in terms of a and b\n')
P_dens_Tc  = P ==   pc_val*kB*T * (1-b*pc_val)^-1 -   a*pc_val^2;
dPdp_Tc     = diff (P_dens,p);
%dPdp_Tc     = 0  == -   kB*T * (pc_val-b)^-2   + 2*a*pc_val^-3;
Tc_eqn_dens = Tc == solve (dPdp_Tc,T);
Tc_val_dens =       solve (dPdp_Tc,T);

dPdV_Tc     = 0  == - N*kB*T * (Vc_val-N*b)^-2 + 2*N^2*a*Vc_val^-3;
Tc_eqn_Vol  = Tc == solve (dPdV_Tc,T);
Tc_val_Vol  =       solve (dPdV_Tc,T);

fprintf(' %s\t %s\n%s\t\t\t\t\t%s\n\n',dPdp_Tc,dPdV_Tc,Tc_eqn_dens,Tc_eqn_Vol)


fprintf('Substitute pc and Tc into P to find Pc in terms of a and b\n')
Pc_eqn_dens = Pc ==   kB*Tc_val_dens * (pc_val-b)^-1   -     a*pc_val^2;
Pc_val_dens =         kB*Tc_val_dens * (pc_val-b)^-1   -     a*pc_val^2;

Pc_eqn_Vol  = Pc == N*kB*Tc_val_dens * (Vc_val-N*b)^-1 - N^2*a*Vc_val^-2;
Pc_val_Vol  =       N*kB*Tc_val_dens * (Vc_val-N*b)^-1 - N^2*a*Vc_val^-2;

fprintf('%s\t%s\n\n',Pc_eqn_dens,Pc_eqn_Vol)


% And now find T, p and P relative to their critical values.
Tr_eqn_dens = Tr == T/Tc_val_dens;
%Tr_val =            T/Tc_val; % Not used
T_eqn_dens  = T  == solve (Tr_eqn_dens,T);
T_val_dens  =       solve (Tr_eqn_dens,T);

Tr_eqn_Vol  = Tr == T/Tc_val_dens;
%Tr_val =            T/Tc_val; % Not used
T_eqn_Vol   = T  == solve (Tr_eqn_Vol,T);
T_val_Vol   =       solve (Tr_eqn_Vol,T);


pr_eqn      = pr == p/pc_val;
%pr_val      =       p/pc_val; % Not used
p_eqn       = p  == solve (pr_eqn,p);
p_val       =       solve (pr_eqn,p);

Vr_eqn      = Vr == V/Vc_val;
%pr_val      =       p/pc_val; % Not used
V_eqn       = V  == solve (Vr_eqn,V);
V_val       =       solve (Vr_eqn,V);


Pr_eqn_dens = Pr == P/Pc_val_dens;
%Pr_val =            P/Pc_val; % Not used
P_eqn_dens  = P  == solve (Pr_eqn_dens,P);
P_val_dens  =       solve (Pr_eqn_dens,P);

Pr_eqn_Vol  = Pr == P/Pc_val_Vol;
%Pr_val =            P/Pc_val; % Not used
P_eqn_Vol   = P  == solve (Pr_eqn_Vol,P);
P_val_Vol   =       solve (Pr_eqn_Vol,P);


Tr_eqn_dens = Tr == T/Tc;
Tr_eqn_Vol  = Tr == T/Tc;

pr_eqn = pr == p/pc;
Vr_eqn = Vr == V/Vc;

Pr_eqn_dens = Pr == P/Pc;
Pr_eqn_Vol  = Pr == P/Pc;

fprintf('%s\t%s\n%s\t%s\n%s\t%s\n\n%s\t%s\n%s\t\t\t\t%s\n%s\t%s\n\n',Tr_eqn_dens,Tr_eqn_Vol,pr_eqn,Vr_eqn,Pr_eqn_dens,Pr_eqn_Vol,T_eqn_dens,T_eqn_Vol,p_eqn,V_eqn,P_eqn_dens,P_eqn_Vol)


P_eqn_dens    = P_val_dens ==   p_val*kB*T_val_dens * (1-b*p_val)^-1    -     a*p_val^2;
%                        P ==   (   N*T*kB   ) / ( V - N*b )    - (N^2*a) / V^2
P_eqn_Vol     = P_val_Vol  ==   N*kB*T_val_Vol * (V_val-N*b)^-1 - N^2*a*V_val^-2;

P_eqn_dens = Pr == solve (P_eqn_dens,Pr);
P_eqn_Vol  = Pr == solve (P_eqn_Vol ,Pr);
P_val_dens =       solve (P_eqn_dens,Pr);
P_val_Vol  =       solve (P_eqn_Vol ,Pr);

fprintf('The Dimensionless Van Der Waals Equation:\n%s\t%s\n\n',P_eqn_dens,P_eqn_Vol)
fprintf('Equivalent to Pr = 8Tr/(3?r-1) - 3?^-2\n')

%% -- Convert to V from p . Not working.

%dens_eqn = p_val == N/V;
%p_val_Vol = solve (dens_eqn,pr);
%
%%pc_eqn_Vol  = pc == N/V;
%%Vc_eqn_dens = Vc == solve (pc_eqn_Vol,V);
%%Vc_val_dens =       solve (pc_eqn_Vol,V);
%
%P_eqn_Vol = P_val_dens == kB*T_val_dens * (p_val_Vol-b)^-1 - a*p_val_Vol^-2;
%P_eqn_Vol = Pr         == solve (P_eqn_Vol,Pr);
%P_val_Vol =               solve (P_eqn_dens,Pr);
%
%%fprintf('%s\n\n',P_eqn_Vol)
%
%%Define VDW Equations in terms of Volume as the script doesn't work so
%%far.
%P_eqn_Vol = Pr == 8*Tr*(3*Vr-1)^-1 - 3*Vr^-2;
%P_val_Vol =      8*Tr*(3*Vr-1)^-1 - 3*Vr^-2;
%
%%fprintf('%s\n\n',P_eqn_Vol)
%
%%-- End convert to V from p.

%% Energy

fprintf('Calculating Energy Now\n')

Vc_val = 3*N*b; V_val  = 3*N*b*Vr; % I haven't got the VDW
%                               equation working just yet in
%                               terms of Volume, so am defining
%                               the values here.

E_eqn_Vol  = E == 3*N*kB*Tc_val_dens/2 - N^2*a/Vc_val;

Ec_eqn_Vol = Ec == solve (E_eqn_Vol,E);
Ec_val_Vol =       solve (E_eqn_Vol,E);

Er_eqn_Vol = Er == E/Ec_val_Vol;
E_eqn_Vol  = E  == solve (Er_eqn_Vol,E);
E_val_Vol  =       solve (Er_eqn_Vol,E);

fprintf('%s\n\n',E_eqn_Vol)


E_eqn_Vol  = E_val_Vol == 3*N*kB*T_val_dens/2 - N^2*a/V_val;
E_eqn_Vol  = Er == solve (E_eqn_Vol,Er);
E_val_Vol  =       solve (E_eqn_Vol,Er);


fprintf('%s\n\n',E_eqn_Vol)

%% Enthalpy

fprintf('Calculating Enthalpy\n')

Hc_eqn_Vol = Hc == Ec_val_Vol + P_val_Vol*Vc_val;
Hc_val_Vol =       Ec_val_Vol + P_val_Vol*Vc_val;

Hr_eqn_Vol = Hr == H/Hc_val_Vol;
H_eqn_Vol  = H  == solve (Hr_eqn_Vol,H);
H_val_Vol  =       solve (Hr_eqn_Vol,H);

fprintf('%s\n\n',H_eqn_Vol)


H_eqn_Vol  = H_val_Vol == E_val_Vol + P_val_Vol*V_val;
H_eqn_Vol  = Hr        == solve (H_eqn_Vol,Hr);
%H_val_Vol  =              solve (H_eqn_Vol,Hr); % Not used

fprintf('%s\n\n',H_eqn_Vol)



%% Gibbs Energy

fprintf('Calculating Gibbs Free Energy\n')

%A_dens = -1 * int(P_val_dens,pr) * p_val/V_val
%A_Vol  = -1 * int(P_val_Vol, Vr);

%G_eqn_Vol = G == A_Vol + P_val_Vol*V_val;

Gc_eqn_Vol = Gc == solve (G_eqn_Vol,G);
Gc_val_Vol =       solve (G_eqn_Vol,G);

Gr_eqn_Vol = Gr == G/Gc_val_Vol
G_eqn_Vol  = G  == solve (Gr_eqn_Vol,G)
G_val_Vol  =       solve (Gr_eqn_Vol,G)

fprintf('%s\n\n',G_eqn_Vol)


%G_eqn_Vol  = G_val_Vol == 3*N*kB*T_val_Vol/2 - N^2*a/V_val;
%G_eqn_Vol  = G_val_Vol == - 3/V_val - (8*T_val_Vol*log(V_val - 1/3))/3 + P_val_Vol*V_val
%G_eqn_Vol  = G_val_Vol == int(P_val_Vol,Vr) + P_val_Vol*V_val
G_eqn_Vol  = G_val_Vol == int(Pc_val_Vol*P_val_Vol,Vr) + Pc_val_Vol*P_val_Vol*V_val
G_eqn_Vol  = Gr == solve (G_eqn_Vol,Gr)
G_val_Vol  =       solve (G_eqn_Vol,Gr)



-Gr*(3/Vr + (8*Tr*log(Vr - 1/3))/3 + (3*N*Vr*b*(8*Tr*Vr^2 - 9*Vr + 3))/(- 3*Vr^3 + Vr^2)) == - 1/(N*Vr*b) - (3*N*Vr*b*(8*Tr*Vr^2 - 9*Vr + 3))/(- 3*Vr^3 + Vr^2) - (64*Tr*a*log(3*N*Vr*b - 1/3))/(81*b*kB)




fprintf('%s\n\n',G_eqn_Vol)

%end