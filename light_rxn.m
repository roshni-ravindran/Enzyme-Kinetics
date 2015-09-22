%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Fernando Contreras
% @email: f2contre@gmail.com
% @project: FIAT LUX
% @institution: University of California, San Diego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system of equations for the light production pathway


function dydt = light_rxn(t,y,P,C)

% t: time
% y:individual variables 
% P: parameters
% C: fixed concentrations

%enzymes 
luxAB = y(1);
frp = y(2);
%complexes
luxAB_FMNH2 = y(3);
luxAB_FMNH2O2 = y(4);
luxAB_FMNH2O2_RCHO = y(5);
luxAB_FMNHOH_RCHO = y(6);
luxAB_FMNH2O2_RCHO_excited = y(7);
luxAB_FMNHOH = y(8);
luxAB_RCHO = y(9);
%products
FMNH2 = y(10);
FMN = y(11);
RCHO = y(12);

%parameters, rate constants
k7  = P(1);
k8  = P(2);
k9  = P(3);
k10  = P(4);
k11  = P(5);
k12  = P(6);
k13  = P(7);
k14  = P(8);
k15  = P(9);
k16 = P(10);
k17 = P(11);
k18 = P(12);
k19 = P(13);
k20 = P(14);
k21 = P(15);

%fixed concentrations
NADPH = C(1);
H = C(2);
O2 = C(3);

dydt = zeros(13,1); %initiate dydt vector

%rate equations
dydt(1) = k10*luxAB_FMNH2 - k9*luxAB*FMNH2 + k12*luxAB_FMNH2O2 + k19*luxAB_FMNHOH + k21*luxAB_RCHO - k20*luxAB*RCHO;%dluxAB(t)/dt
dydt(2) = k7*NADPH*H*FMN*frp-k7*NADPH*H*FMN*frp;%dfrp(t)/dt %constant concentration, slow degradation
dydt(3) = k9*luxAB*FMNH2 - k10*luxAB_FMNH2 - k11*luxAB_FMNH2*O2;%dluxAB-FMNH2(t)/dt
dydt(4) = k11*luxAB_FMNH2*O2 - k12*luxAB_FMNH2O2 + k14*luxAB_FMNH2O2_RCHO - k13*RCHO*luxAB_FMNH2O2;%dluxAB_FMNH2O2(t)/dt
dydt(5) = k13*RCHO*luxAB_FMNH2O2 - k15*luxAB_FMNH2O2_RCHO - k14*luxAB_FMNH2O2_RCHO;%dluxAB-FMNH2O2-RCOH(t)/dt
dydt(6) = k17*RCHO*luxAB_FMNHOH - k18*luxAB_FMNHOH_RCHO;%dluxAB-FMNHOH-RCOH(t)/dt
dydt(7) = k15*luxAB_FMNH2O2_RCHO - k16*luxAB_FMNH2O2_RCHO_excited;%dluxAB-FMNH2O2-RCOH*(t)/dt
dydt(8) = k16*luxAB_FMNH2O2_RCHO_excited + k18*luxAB_FMNHOH_RCHO - k17*RCHO*luxAB_FMNHOH - k19*luxAB_FMNHOH;%dluxAB-FMNHOH(t)/dt
dydt(9) = k20*luxAB*RCHO - k21*luxAB_RCHO;%dluxAB-RCOH(t)/dt
dydt(10) = k7*NADPH*H*FMN*frp - k8*O2*FMNH2 - k9*luxAB*FMNH2 + k10*luxAB_FMNH2;%dFMNH2(t)/dt
dydt(11) = k8*O2*FMNH2 - k7*NADPH*H*FMN*frp + k12*luxAB_FMNH2O2 + k19*luxAB_FMNHOH;%dFMN(t)/dt
dydt(12) = k14*luxAB_FMNH2O2_RCHO - k13*RCHO*luxAB_FMNH2O2 + k18*luxAB_FMNHOH_RCHO - k17*RCHO*luxAB_FMNHOH + k21*luxAB_RCHO - k20*luxAB*RCHO;%dRCOH(t)/dt
dydt(13) = k16*luxAB_FMNH2O2_RCHO_excited;%dhv(t)/dt, dRCOOH(t)/dt