%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Fernando Contreras
% @email: f2contre@gmail.com
% @project: FIAT LUX
% @institution: University of California, San Diego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Can be utilized to effectively study single-burst assays in vitro and the 
%sensitivity of certain parameters within the light production pathway

clc; 
clear;
tic
micron = 10^(-6);

%enzymes
luxAB = 10*micron;
frp = 10*micron;
%enzyme-substrate complexes
luxAB_FMNH2 = 0;
luxAB_FMNH2O2 = 0;
luxAB_FMNH2O2_RCOH = 0;
luxAB_FMNHOH_RCOH = 0;
luxAB_FMNH2O2_RCOH_excited = 0;
luxAB_FMNHOH = 0;
luxAB_RCOH = 0;
%substrates and precursors
FMNH2 = 100*micron; 
FMN = 0*micron;
RCHO = 100*micron;
hv = 0;

%initial conditions
y0 = [luxAB,frp,luxAB_FMNH2,luxAB_FMNH2O2,luxAB_FMNH2O2_RCOH,luxAB_FMNHOH_RCOH,...
    luxAB_FMNH2O2_RCOH_excited,luxAB_FMNHOH,luxAB_RCOH,FMNH2,FMN,RCHO,hv];

%environmental conditions
NADPH = 560*micron;
H = 300*micron;
O2 = 550*micron;

%fixed concentrations
C = [NADPH,H,O2];

%rate constants
k7  = 21.2;%s^(-1)
k8  = 10;%s^(-1)
k9  = 6*10^5;%M^(-1)s^(-1)
k10  = 4.6;%s^(-1)
k11  = 2.4*10^6;%M^(-1)s^(-1)
k12  = 0.10;%s^(-1)
k13  = 1.2*10^5;%M^(-1)s^(-1)
k14  = 0.1;%s^(-1)
k15  = 9.5;%s^(-1)
k16 = 0.5;%s^(-1)
k17 = 3*10^3;%M^(-1)s^(-1)
k18 = 0.06;%s^(-1)
k19 = 1.9*10^-3;%s^(-1)
k20 = 1*10^5;%M^(-1)s^(-1)
k21 = 40;%s^(-1) 

%parameters
P = [k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21];

% tspan = 0:0.05:100; %simulation time
% [t,y] = ode23(@light_rxn,tspan,y0,[],P,C);
% light = y(:,end);

% figure;
% pltAxis = gca;
% hold(pltAxis,'on')
% xlabel('Time (s)','Fontsize',12)
% ylabel('hv (au)' ,'Fontsize',12)
% set(pltAxis,'Fontsize',12)
% 
% plot(pltAxis,t,light./micron,'Linewidth',3,'color','b')
% toc

%% [RCHO] vs hv_ss

clc; 
tic

%initiate figure for individual traces
s = figure;
pltAxis = gca; 
set(s,'name','RCHO vs hv','numbertitle','off')
hold(pltAxis,'on')
xlabel('Time (s)','Fontsize',15)
ylabel('Light (a.u.)','Fontsize',15)
set(pltAxis,'Fontsize',15)
% set(pltAxis,'xlim',[0 15])

%parameter space: two orders of magnitude above and below the nominal parameter value
order = 2;
pSpace = RCHO*(10.^(-order:order)); % [param](uM)

cmap = hsv(numel(pSpace)); %traces colors
tspan = 0:0.05:20; %simulation time

%initiate vectors for parameter vs steady light figure
param_rcho = pSpace./micron;
ss_rcho = zeros(numel(tspan),numel(pSpace));

%legend
p = zeros(1,numel(pSpace));

for i = 1:numel(pSpace)
    y0(12) = pSpace(i);
    
    [t,y] = ode23(@light_rxn,tspan,y0,[],P,C);
    light = y(:,end)./micron;
    p(i) = plot(pltAxis,t,light,'Linewidth',3,'color',cmap(i,:));
    
    ss_rcho(:,i) = light;
end
y0(12) = RCHO;

legend(p,'[RCHO] = 1*10^-^6 M','[RCHO] = 1*10^-^5 M','[RCHO] = 1*10^-^4 M',...
    '[RCHO] = 1*10^-^3 M','[RCHO] = 1*10^-^2 M','Location','best');

figure;
x = log10(param_rcho);
plot(x,ss_rcho(end,:),'*-','Linewidth',2,'color','k')
xlabel('log_1_0([RCHO]) (\muM)', 'Fontsize',15)
ylabel('Steady State Light (a.u.)', 'Fontsize',15)
set(gca,'Fontsize',15)
set(gca,'xlim',[x(1) x(end)])
toc


%% [FMNH2] vs hv_ss

clc; 
tic

%initiate figure for individual traces
s = figure;
pltAxis = gca; 
set(s,'name','FMNH2 vs hv','numbertitle','off')
hold(pltAxis,'on')
xlabel('Time (s)','Fontsize',15)
ylabel('Light (a.u.)','Fontsize',15)
set(pltAxis,'Fontsize',15)
% set(pltAxis,'xlim',[0 15])

%parameter space: two orders of magnitude above and below the nominal parameter value
order = 2;
pSpace = FMNH2*(10.^(-order:order)); % [param](uM)

cmap = hsv(numel(pSpace)); %traces colors
tspan = 0:0.05:20; %simulation time

%initiate vectors for parameter vs steady light figure
param_fmnh2 = pSpace./micron;
ss_fmnh2 = zeros(numel(tspan),numel(pSpace));

%legend
p = zeros(1,numel(pSpace));

for i = 1:numel(pSpace)
    y0(10) = pSpace(i);
    
    [t,y] = ode23(@light_rxn,tspan,y0,[],P,C);
    light = y(:,end)./micron;
    p(i) = plot(pltAxis,t,light,'Linewidth',3,'color',cmap(i,:));
    
    ss_fmnh2(:,i) = light;
end
y0(10) = FMNH2;

legend(p,'[FMNH_2] = 1*10^-^6 M','[FMNH_2] = 1*10^-^5 M','[FMNH_2] = 1*10^-^4 M',...
    '[FMNH_2] = 1*10^-^3 M','[FMNH_2] = 1*10^-^2 M','Location','best');

figure;
x = log10(param_fmnh2);
plot(x,ss_fmnh2(end,:),'*-','Linewidth',2,'color','k')
xlabel('log_1_0([FMNH_2]) (\muM)', 'Fontsize',15)
ylabel('Steady State Light (a.u.)', 'Fontsize',15)
set(gca,'Fontsize',15)
set(gca,'xlim',[x(1) x(end)])

toc

%% compare FMNH2 and RCHO steady state light production 
% figure;
% x = log10(param_fmnh2);
% hold on
% p1 = plot(x,ss_fmnh2(end,:),'*-','Linewidth',2,'color','r');
% p2 = plot(x,ss_rcho(end,:),'*-','Linewidth',2,'color','b');
% xlabel('log_1_0([ ]) (M)', 'Fontsize',12)
% ylabel('hv_s_s (AU)', 'Fontsize',12)
% set(gca,'Fontsize',12)
% set(gca,'xlim',[x(1) x(end)])
% 
% legend([p1,p2],'FMNH_2','RCHO','Location','best')


%% k19 vs hv_ss
% 
clc; 
tic

%initiate figure for individual traces
s = figure;
pltAxis = gca; 
set(s,'name','k19 vs hv','numbertitle','off')
hold(pltAxis,'on')
xlabel('Time (s)','Fontsize',15)
ylabel('Light (a.u.)','Fontsize',15)
set(pltAxis,'Fontsize',15)
% set(pltAxis,'xlim',[0 15])

%parameter space: two orders of magnitude above and below the nominal parameter value
order = 3;
pSpace = 1*(10.^(-order:order)); % [param](uM)

cmap = hsv(numel(pSpace)); %traces colors
tspan = 0:0.05:60; %simulation time

%initiate vectors for parameter vs steady light figure
param = pSpace./micron;
ss = zeros(numel(tspan),numel(pSpace));

%legend
p = zeros(1,numel(pSpace));

for i = 1:numel(pSpace)
    P(13) = pSpace(i);
    
    [t,y] = ode23(@light_rxn,tspan,y0,[],P,C);
    light = y(:,end)./micron;
    p(i) = plot(pltAxis,t,light,'Linewidth',3,'color',cmap(i,:));
    
    ss(:,i) = light;
end
P(13) = k19;

legend(p,'k_1_9 = 1*10^-^3 sec^-^1','k_1_9 = 1*10^-^2 sec^-^1','k_1_9 = 1*10^-^1 sec^-^1',...
    'k_1_9 = 1*10^0 sec^-^1','k_1_9 = 1*10^1 sec^-^1','k_1_9 = 1*10^2 sec^-^1','k_1_9 = 1*10^3 sec^-^1',...
    'Location','best');

figure;
x = log10(param);
plot(x,ss(end,:),'*-','Linewidth',2,'color','k')
xlabel('log_1_0(k_1_9) (\musec^-^1)', 'Fontsize',15)
ylabel('Steady State Light (a.u.)', 'Fontsize',15)
set(gca,'Fontsize',15)
set(gca,'xlim',[x(1) x(end)])
toc

%% [O2] vs hv_ss

% clc; 
% tic
% 
% %initiate figure for individual traces
% s = figure;
% pltAxis = gca; 
% set(s,'name','O2 vs hv','numbertitle','off')
% hold(pltAxis,'on')
% xlabel('Time (s)','Fontsize',15)
% ylabel('hv (au)','Fontsize',15)
% set(pltAxis,'Fontsize',15)
% set(pltAxis,'xlim',[0 8])
% 
% %parameter space: two orders of magnitude above and below the nominal parameter value
% order = 2;
% pSpace = O2*(10.^(-order:order)); % [param](uM)
% 
% cmap = hsv(numel(pSpace)); %traces colors
% tspan = 0:0.05:10; %simulation time
% 
% %initiate vectors for parameter vs steady light figure
% param_fmh2 = pSpace;
% ss_fmnh2 = zeros(numel(tspan),numel(pSpace));
% 
% %legend
% p = zeros(1,numel(pSpace));
% 
% for i = 1:numel(pSpace)
%     C(3) = pSpace(i);
%     
%     [t,y] = ode23(@light_rxn,tspan,y0,[],P,C);
%     light = y(:,end)./micron;
%     p(i) = plot(pltAxis,t,light,'Linewidth',3,'color',cmap(i,:));
%     
%     ss_fmnh2(:,i) = light;
% end
% C(3) = O2;
% 
% legend(p,'[O_2] = 5.5*10^-^6','[O_2] = 5.5*10^-^5','[O_2] = 5.5*10^-^4',...
%     '[O_2] = 5.5*10^-^3','[O_2] = 5.5*10^-^2','Location','best');
% 
% figure;
% x = log10(param_fmh2./micron);
% plot(x,round(ss_fmnh2(end,:)),'*-','Linewidth',2,'color','k')
% xlabel('log_1_0([O2]) (M)', 'Fontsize',12)
% ylabel('hv_s_s (AU)', 'Fontsize',12)
% set(gca,'Fontsize',12)
% set(gca,'xlim',[x(1) x(end)])
% toc