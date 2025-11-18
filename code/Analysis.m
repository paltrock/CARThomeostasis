%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analyze Kimmel model with parameters from fit
%
% https://doi.org/10.1098/rspb.2021.0229
%
% dNdt = -rN*N*ln(T/kN)
% dCdt = -rc(T)*C*ln(T/kC)
% dBdt = rB*B-gB*B*C/(C+kB)
%
% with 
% T = N+C
% rC(T)=pC+b(T-kN)^2/(aT^2+(T-KN)^2)
%
% ------------------------------------------------------------------------

%% Data

Data=readtable("/Users/alvaro/Desktop/Full_dataset.xlsx");
Data([4,15,17,18],:)=[]; % Remove KTE 6,14,20,21 
ALC=Data{:,18:25}*1000*5e6; % 5L=5e6uL (ref?)
CART=Data{:,26:36}*5e6; % 5L=5e6uL (ref?)
Tumor_initial=Data{:,14}*1e8; % 1mL=1cm3=1e8 cells (del monte 2009, wikipedia)

Tumor_response=readmatrix("/Users/alvaro/Desktop/Tumor_response.xlsx");
Tumor_response=Tumor_response(:,2:4);
Tumor_response2=double(Tumor_response~=0);
IDs=Data{:,1}; % Patients' IDs

Parameters=readmatrix("/Users/alvaro/Desktop/Full_parameters.xlsx");
Parameters=Parameters(:,2:end);
Parameters(:,[1:2,7:8])=Parameters(:,[1:2,7:8])*5e6; 

tint=0.1;
tspan=0:tint:180;

Times_PET=[30 90 180];

%% Correlations

figure('Position',[10 10 800 800])
[R,Pval]=corrplot(Parameters,VarNames={'kN','kC','rN','pC','a','b','N0','C0','c','rB','gB','kB'});
set(findall(gcf, 'type', 'axes'), 'XTickLabel', [])
set(findall(gcf, 'type', 'axes'), 'YTickLabel', [])
title([])
%exportgraphics(gcf,'Correlations.png','Resolution',300)

figure
imagesc(R)
colorbar
xticks(1:12)
xticklabels({'kN','kC','rN','pC','a','b','N0','C0','c','rB','gB','kB'})
yticks(1:12)
yticklabels({'kN','kC','rN','pC','a','b','N0','C0','c','rB','gB','kB'})
%exportgraphics(gcf,'Correlations_2.png','Resolution',300)


figure
imagesc(Pval)
colorbar
xticks(1:12)
xticklabels({'kN','kC','rN','pC','a','b','N0','C0','c','rB','gB','kB'})
yticks(1:12)
yticklabels({'kN','kC','rN','pC','a','b','N0','C0','c','rB','gB','kB'})
%exportgraphics(gcf,'Pval.png','Resolution',300)



%% Compute AUC, Cmax and time to recurrence

% Initialize

AUC=zeros(length(IDs),60);
Cmax=zeros(length(IDs),60);
AUC_m=zeros(length(IDs),1);
AUC_p=zeros(length(IDs),1);
Cmax_m=zeros(length(IDs),1);
Cmax_p=zeros(length(IDs),1);
ttr_m=zeros(length(IDs),1);
ttr_p=zeros(length(IDs),1);

% Solve for each patient

for patient=1:length(IDs)

% Resolution ---------------------

param=[Parameters(patient,[1:6,9:12])];

N0=Parameters(patient,7);
C0=Parameters(patient,8);
B0=Tumor_initial(patient);
n0=[N0,C0,B0];

Opt    = odeset('Events', @myEvent);
[t,n] = ode45(@(t,n) model(t,n,param),tspan,n0',Opt);

N=n(:,1);
C=n(:,2);
B=n(:,3);

% ----------------------------------

% AUC & Cmax - complete

for day=1:60

idx=day/tint;
AUC(patient,day)=trapz(C(1:idx));
Cmax(patient,day)=max(C);

end

% -------------

if isnan(CART(patient,2))==0 % Distinguish patients who lack data on day 7

    % AUC & Cmax - measured

    AUC_m(patient)=trapz(CART(patient,[1 2 3 4]));
    Cmax_m(patient)=max(CART(patient,[1 2 3 4]));

    % AUC & Cmax - predicted

    AUC_p(patient)=trapz(C([0 7 14 28]/tint +1));
    Cmax_p(patient)=max(C([0 7 14 28]/tint +1));

else

    AUC_m(patient)=trapz(CART(patient,[1 3 4]));
    AUC_p(patient)=trapz(C([0 14 28]/tint +1));

    Cmax_m(patient)=max(CART(patient,[1 2 3 4]));
    Cmax_p(patient)=max(C([0 7 14 28]/tint +1));

end

% Time to recurrence - measured

Response_clinic=Tumor_response(patient,:);
idx=find(Response_clinic>1);
if isempty(idx)==0

    ttr_m(patient)=Times_PET(idx(1));

else
    
    ttr_m(patient)=0;

end



% Time to recurrence - predicted

[Bmin idxmin]=min(B);
relapse=find(B(idxmin:end)>min(B) & B(idxmin:end)>1e8);

if isempty(relapse)==0
    ttr_p(patient)=(idxmin+relapse(1))*tint;

else
    ttr_p(patient)=0;
end

end





%% Boxplots

figure('Position', [10 10 800 350])
bar(log10([AUC_m AUC_p AUC(:,28)]))
xticklabels([])
xlabel('Patients')
ylabel('log_{10}(AUC-28)')
set(gca,'Fontsize',20)
legend('Measured','Predicted','Complete','Location','northoutside','Orientation','horizontal')
%exportgraphics(gcf,['AUC-28_Comparison.png'],'Resolution',300)

figure('Position', [10 10 800 350])
bar(log10([Cmax_m Cmax_p Cmax(:,28)]))
xticklabels([])
xlabel('Patients')
ylabel('log_{10}(C_{max})')
set(gca,'Fontsize',20)
legend('Measured','Predicted','Complete','Location','northoutside','Orientation','horizontal')
%exportgraphics(gcf,['Cmax_Comparison.png'],'Resolution',300)


day=28;
% Kolmogorov-Smirnov test
[h,p] = kstest2(AUC(logical(Tumor_response2(:,3)),day),AUC(~logical(Tumor_response2(:,3)),day));
%[h,p] = kstest2(AUC_a(logical(Tumor_response2(:,3))),AUC_a(~logical(Tumor_response2(:,3))));

figure
b1=boxplot(AUC(:,day),Tumor_response2(:,3),'labels',{'Response','No response'});
%b1=boxplot(AUC_a,Tumor_response2(:,3),'labels',{'Response','No response'});
%xlabel('Response at day 180')
ylabel(['AUC-',num2str(day)])
annotation('textbox',[.6 .5 .3 .3],'String',strcat("p-val= ",num2str(p)),'FitBoxToText','on','fontsize',14)
set(b1,'linewidth',2)
set(gca,'FontSize',20,'linewidth',2)
%exportgraphics(gcf,['AUC-',num2str(day),'.png'],'Resolution',300)


for day=1:60
[h,p(day)] = kstest2(AUC(logical(Tumor_response2(:,3)),day),AUC(~logical(Tumor_response2(:,3)),day));
end

figure
plot(1:60,p,'linewidth',3)
xlabel('AUC day')
ylabel('p-value')
hold on
yline(0.05,'r--','linewidth',3)
set(gca,'FontSize',20)
%exportgraphics(gcf,'pvalAUC.png','Resolution',300)


%[h,p] = kstest2(Cmax(logical(Tumor_response2(:,3)),day),Cmax(~logical(Tumor_response2(:,3)),day));
[h,p] = kstest2(Cmax_p(logical(Tumor_response2(:,3))),Cmax_p(~logical(Tumor_response2(:,3))));

figure
%b1=boxplot(Cmax(:,day),Tumor_response2(:,3),'labels',{'Response','No response'});
b1=boxplot(Cmax_p,Tumor_response2(:,3),'labels',{'Response','No response'});
%xlabel('Response at day 180')
ylabel('CAR peak')
annotation('textbox',[.6 .5 .3 .3],'String',strcat("p-val= ",num2str(p)),'FitBoxToText','on','fontsize',14)
set(b1,'linewidth',2)
set(gca,'FontSize',20,'linewidth',2)
%exportgraphics(gcf,'Cmax_a.png','Resolution',300)



%% Plot ttr

% Find day of recurrence
idx=find(ttr_m~=0);
% Set dimensions for scatter plot
sz=2000*Tumor_initial(idx)/max(Tumor_initial);


figure
scatter(ttr_m(idx),AUC(idx),sz,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980])
xlabel('Time to recurrence (measured)')
ylabel('AUC-30')
xlim([0,200])
set(gca,'Fontsize',20)
%exportgraphics(gcf,'AUC-ttrm.png','Resolution',300)


figure
scatter(ttr_p(idx),AUC(idx),sz,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','flat')
xlabel('Time to recurrence (predicted)')
ylabel('AUC-30')
xlim([0,200])
set(gca,'Fontsize',20)
%exportgraphics(gcf,'AUC-ttrp.png','Resolution',300)

figure
scatter(ttr_p(idx),AUC(idx),sz,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','flat')
hold on
scatter(ttr_m(idx),AUC(idx),sz,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','flat')
xlabel('Time to recurrence (combined)')
ylabel('AUC-30')
xlim([0,200])
legend('Predicted','Measured')
set(gca,'Fontsize',20)
%exportgraphics(gcf,'AUC-ttr.png','Resolution',300)


%% Functions

% Model

function dndt = model(t,n,param)

dndt=zeros(3,1);

kN=param(1);    
kC=param(2);   
rN=param(3);      
pC=param(4);   
a=param(5);    
b=param(6); 
c=param(7);
rB=param(8);     
gB=param(9);      
kB=param(10);   

T=n(1)+n(2);
rc=pC+b*((T-kN)^2)/(a*T^2+(T-kN)^2);
dndt(1)=-rN*n(1)*log((n(1)+c*n(2))/kN);
dndt(2)=-rc*n(2)*log(T/kC);
dndt(3)=rB*n(3)-gB*n(3)*n(2)/(n(2)+kB);

if n(1)<1
    dndt(1)=0;
end

if n(2)<1
    dndt(2)=0;
end

if n(3)<1
    dndt(3)=0;
end

end

function [value, isterminal, direction] = myEvent(T, Y)
value      = (Y(3) > 1e12);
isterminal = 1;   % Stop the integration
direction  = 0;
end

% Compute response at days 30,90,180

function R=response(B,t)

tint=t(2)-t(1);
R=zeros(1,3);
T=length(t);

        idx30=30/tint;
        if idx30<T
 
            if B(idx30)<1e8
                R(1)=0;
            elseif B(idx30)<0.9*B(1)
                R(1)=1;
            elseif B(idx30)<1.1*B(1)
                R(1)=2;
            else
                R(1)=3;
            end

        else
            R(1)=4;
        end

    

        idx90=90/tint;
        if idx90<T

            if B(idx90)<1e8
                R(2)=0;
            elseif B(idx90)<0.95*B(1)
                R(2)=1;
            elseif B(idx90)<1.05*B(1)
                R(2)=2;
            else
                R(2)=3;
            end

        else
            R(2)=4;
        end



        idx180=180/tint;
        if idx180<T

            if B(idx180)<1e8
                R(3)=0;
            elseif B(idx180)<0.95*B(1)
                R(3)=1;
            elseif B(idx180)<1.05*B(1)
                R(3)=2;
            else
                R(3)=3;
            end

        else
            R(3)=4;
        end
        


end