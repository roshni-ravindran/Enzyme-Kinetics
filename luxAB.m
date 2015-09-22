%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Fernando Contreras
% @email: f2contre@gmail.com
% @project: FIAT LUX
% @institution: University of California, San Diego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%script used to simulate ligth emission over time

function dydt = luxAB(t,y,C)

FMNH2 = y(1);
RCHO = y(2);

%parameters
Vmax_luxAB = 71.58;
Ki_FMNH2 = 0.62;
Km_O2 = 81.5;
Km_FMNH2 = 0.22;
Km_RCHO = 72.2;

%fixed concentration 
O2 = C(1);

v_luxAB = (Vmax_luxAB*(FMNH2)*(O2)*(RCHO))/(Ki_FMNH2*Km_O2*(RCHO) + Km_FMNH2*(O2).*(RCHO) + Km_O2*(FMNH2)*(RCHO) + Km_RCHO*(FMNH2)*(O2) + (FMNH2)*(O2)*(RCHO));

dydt = zeros(3,1);
dydt(1) = -v_luxAB;
dydt(2) = -v_luxAB;
dydt(3) = v_luxAB;

