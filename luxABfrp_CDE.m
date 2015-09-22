%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Fernando Contreras
% @email: f2contre@gmail.com
% @project: FIAT LUX
% @institution: University of California, San Diego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%script used to simulate the behavior of our control 

function dydt = luxABfrp_CDE(t,y,C,P)

FMNH2 = y(1);
RCHO = y(2);
FMN = y(3);
RCOOH = y(4);
RCOACP = y(5);

%frp
Ki_FMN = 1.0;
Km_NADPH = 49.5;
Km_FMN = 0.72;

%luxAB
Ki_FMNH2 = 0.62;
Km_O2 = 81.5;
Km_FMNH2 = 0.22;
Km_RCHO = 72.2;

%luxEC
r44 = 0.04;
k62 = 95.3;
k61 = 90.9;
k63 = 24.35;
k64 = 76.5;

%luxD
Ki_RCOACP = P(1);
Km_RCOACP = 0.37;
Km_H2O = P(2);

%Vmax
Vmax_frp = 51.8;
Vmax_luxAB = 71.58;
Vmax_luxEC = 198.93;
Vmax_luxD = 45.98;

%fixed concentrations (uM)
NADPH = 560;
O2 = 550;
ATP = 1310;
H2O = 500;

%enzyme concentrations
luxAB = C(1);
frp = C(2);
luxEC = C(3);
luxD = C(4);

tau = 0.067;

v_frp = (Vmax_frp*(FMN)*(NADPH))./(Ki_FMN*Km_NADPH + Km_NADPH*(FMN) + Km_FMN*(NADPH) + (FMN)*(NADPH));
v_luxAB = (Vmax_luxAB*(FMNH2)*(O2)*(RCHO))./(Ki_FMNH2*Km_O2*(RCHO) + Km_FMNH2*(O2)*(RCHO) + Km_O2*(FMNH2)*(RCHO) + Km_RCHO*(FMNH2)*(O2) + (FMNH2)*(O2)*(RCHO));
v_luxEC = (Vmax_luxEC*(RCOOH)*(ATP)*(NADPH))./(r44*k62*(NADPH) + k61*(ATP)*(NADPH) + k62*(RCOOH)*(NADPH) + 2*k63*(RCOOH)*(ATP)*(NADPH)+ 2*k64*(RCOOH)*(ATP) + (RCOOH)*(ATP)*(NADPH));
v_luxD = (Vmax_luxD*(RCOACP)*(H2O))./(Ki_RCOACP*Km_H2O + Km_H2O*(RCOACP) + Km_RCOACP*(H2O) + (RCOACP)*(H2O));

P0 = tau*(RCOOH + RCHO + RCOACP);

dydt = zeros(5,1);

dydt(1) = v_frp*frp - v_luxAB*luxAB; %delta[FMNH2]
dydt(2) = P0 + v_luxEC*luxEC - v_luxAB*luxAB - tau*RCHO; %delta[RCHO]
dydt(3) = v_luxAB*luxAB - v_frp*frp; %delta[FMN]
dydt(4) = -v_luxEC*luxEC + v_luxAB*luxAB + v_luxD*luxD - tau*RCOOH; %delta[RCOOH]
dydt(5) = -v_luxD*luxD; %delta[RCOACP]