%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Fernando Contreras
% @email: f2contre@gmail.com
% @project: FIAT LUX
% @institution: University of California, San Diego
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%script used to simulate and analyze the different genetic permutations

%% Analysis of Genetic Permutations

clear; 
clc;

%initial []s (uM)
fmnh2 = 78;
rcho = 0;
fmn = 10;
rcooh = 80;
rcoacp = 150;

y0 = [fmnh2,rcho,fmn,rcooh,rcoacp];

%enzyme []s (uM)
luxab = 10*0.3;
frp = 0.5; 
luxEC = 0.3;
luxD = 0.3;

C = [luxab,frp,luxEC,luxD];

%unknown parameters
Ki_RCOACP = 1;
Km_H2O = 80;

P = [Ki_RCOACP,Km_H2O];

%simulation time
tspan = 0:0.1:200; 

%numerical solutions
%luxCDE
[t_cde,y_cde] = ode23(@luxABfrp_CDE,tspan,y0,[],C,P);
fmnh2_cde = y_cde(:,1);
fmn_cde = y_cde(:,3);
rcho_cde = y_cde(:,2);
rcooh_cde = y_cde(:,4);

%luxCDED
[t_cded,y_cded] = ode23(@luxABfrp_CDED,tspan,y0,[],C,P);
fmnh2_cded = y_cded(:,1);
fmn_cded = y_cded(:,3);
rcho_cded = y_cded(:,2);
rcooh_cded = y_cded(:,4);

%luxCDEC
[t_cdec,y_cdec] = ode23(@luxABfrp_CDE,tspan,y0,[],C,P);
%luxCDEE
[t_cdee,y_cdee] = ode23(@luxABfrp_CDE,tspan,y0,[],C,P);

%luxCDECE
[t_cdece,y_cdece] = ode23(@luxABfrp_CDECE,tspan,y0,[],C,P);
fmnh2_cdece = y_cdece(:,1);
fmn_cdece = y_cdece(:,3);
rcho_cdece = y_cdece(:,2);
rcooh_cdece = y_cdece(:,4);

%rcho
figure;
plt1 = gca;
set(plt1,'xlim',[0, 40]);
hold(plt1,'on')
p1 = plot(plt1,t_cde,rcho_cde,'color','r','Linewidth',3);
p2 = plot(plt1,t_cde,rcho_cde,'color','r','Linewidth',3);
p3 = plot(plt1,t_cde,rcho_cde,'color','r','Linewidth',3);
p4 = plot(plt1,t_cded,rcho_cded,'color','b','Linewidth',3);
p5 = plot(plt1,t_cdece,rcho_cdece,'color','g','Linewidth',3);
legend([p1,p2,p3,p4,p5],'luxCDE','luxCDE-C','luxCDE-E','luxCDE-D','luxCDE-CE','Location','best')
ylabel('[RCHO] (\muM)','Fontsize',15)
xlabel('Time (s)','Fontsize',15)
set(plt1,'Fontsize',15)

%rcooh
figure;
plt2 = gca;
set(plt2,'xlim',[0 40]);
hold(plt2,'on')
p1 = plot(plt2,t_cde,rcooh_cde,'color','r','Linewidth',3);
p2 = plot(plt2,t_cde,rcooh_cde,'color','r','Linewidth',3);
p3 = plot(plt2,t_cde,rcooh_cde,'color','r','Linewidth',3);
p4 = plot(plt2,t_cded,rcooh_cded,'color','b','Linewidth',3);
p5 = plot(plt2,t_cdece,rcooh_cdece,'color','g','Linewidth',3);
legend([p1,p2,p3,p4,p5],'luxCDE','luxCDE-C','luxCDE-E','luxCDE-D','luxCDE-CE','Location','best')
ylabel('[RCOOH] (\muM)','Fontsize',15)
xlabel('Time (s)','Fontsize',15)
set(plt2,'Fontsize',15)

% %fmnh2
% figure;
% plt3 = gca;
% set(plt3,'xlim',[0 40]);
% hold(plt3,'on')
% p1 = plot(plt3,t_cde,fmnh2_cde,'color','r','Linewidth',3);
% p2 = plot(plt3,t_cded,fmnh2_cded,'color','b','Linewidth',3);
% p3 = plot(plt3,t_cdece,fmnh2_cdece,'color','g','Linewidth',3);
% legend([p1,p2,p3],'luxCDE','luxCDED','luxCDECE')
% ylabel('[FMNH_2]')
% 
% %fmn
% figure;
% plt4 = gca;
% set(plt4,'xlim',[0 40]);
% hold(plt4,'on')
% p1 = plot(plt4,t_cde,fmn_cde,'color','r','Linewidth',3);
% p2 = plot(plt4,t_cded,fmn_cded,'color','b','Linewidth',3);
% p3 = plot(plt4,t_cdece,fmn_cdece,'color','g','Linewidth',3);
% legend([p1,p2,p3],'luxCDE','luxCDED','luxCDECE')
% ylabel('[FMN]')

%% hv production 
clc;

hv = 0;
tspan = 0:0.1:15;
o2 = 550;

%luxCDE
fmnh2CDE = max(fmnh2_cde);
rchoCDE = max(rcho_cde);
[t,yCDE] = ode23(@luxAB,tspan,[fmnh2CDE,rchoCDE,hv],[],o2);
hv_cde= yCDE(:,end);

%luxCDED
fmnh2CDED = max(fmnh2_cded);
rchoCDED = max(rcho_cded);
[~,yCDED] = ode23(@luxAB,tspan,[fmnh2CDED,rchoCDED,hv],[],o2);
hv_cded = yCDED(:,end);

%luxCDECE
fmnh2CDECE = max(fmnh2_cdece);
rchoCDECE = max(rcho_cdece);
[~,yCDECE] = ode23(@luxAB,tspan,[fmnh2CDECE,rchoCDECE,hv],[],o2);
hv_cdece = yCDECE(:,end);

figure;
plt7 = gca;
hold(plt7,'on')
ylabel('Light (au)','Fontsize',15)
xlabel('Time (s)','Fontsize',15)
set(plt7,'Fontsize',15)
p1 = plot(plt7,t,hv_cde,'color','r','Linewidth',3);
p2 = plot(plt7,t,hv_cde,'color','r','Linewidth',3);
p3 = plot(plt7,t,hv_cde,'color','r','Linewidth',3);
p4 = plot(plt7,t,hv_cded,'color','b','Linewidth',3);
p5 = plot(plt7,t,hv_cdece,'color','g','Linewidth',3);
legend([p1,p2,p3,p4,p5],'luxCDE','luxCDE-C','luxCDE-E','luxCDE-D','luxCDE-CE','Location','best')

%% barplot
figure;

CDEss = hv_cde(end);
CDEDss = hv_cded(end);
CDECEss = hv_cdece(end);
CDECss = CDEss;
CDEEss = CDEss;

vals = [CDEss,CDEDss,CDECEss,CDECss,CDEEss];

bar(vals);
ylabel('Steady State Light (au)','Fontsize',15)
set(gca,'XTick',[1, 2, 3, 4, 5])
permutations = {'luxCDE','luxCDE-D','luxCDE-CE','luxCDE-C','luxCDE-E'};
set(gca,'XTickLabel',permutations)
set(gca,'Fontsize',15)


