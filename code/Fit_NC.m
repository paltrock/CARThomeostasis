              % Fit Moffit patient data

Data=readtable("/Users/alvaro/Desktop/Full_dataset.xlsx");
ALC=Data{:,18:25}*1000;
CART=Data{:,26:36};
IDs=Data{:,1};

%% Parse data 

% Time intervals
t_CART=[0,7,14,28,90,180,270,360,450,540,720];
t_ALC=[0,1,5,7,14,28,90,180];
t_common=[0,7,14,28,90,180];

% Data at common times
CART_sub=CART(:,1:6);
ALC_sub=ALC(:,[1,4,5,6,7,8]);
TCELL=ALC_sub-CART_sub;

% Compute median
ALC_median=median(ALC_sub,'omitnan');
TCELL_median=median(TCELL,'omitnan');
CART_median=median(CART_sub,'omitnan');
TCELL_sd=std(TCELL,'omitnan');
CART_sd=std(CART_sub,'omitnan');

Data_median=[TCELL_median',CART_median'];
Data_sd=[TCELL_sd',CART_sd'];

%% Plot median data

figure('Position', [10 10 800 350])
yyaxis left
errorbar(t_common,TCELL_median,TCELL_sd/2,'*-','linewidth',3);
ylabel('ALC (Cells/uL)')
xline(0,'--','color','red','LineWidth',2)
yyaxis right 
errorbar(t_common,CART_median,CART_sd/2,'*-','linewidth',3);
ylabel('CAR-T (Cells/uL)')
xlabel('Time (month)')
xlim([0,180])
xticks([0,14,30,60,90,180])
set(gca,'fontsize',20)

%% Fit median data

% Fmincon options
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'Display','iter','StepTolerance',1e-10,'MaxFunctionEvaluations',9000);

% Define upper and lower bounds
lb=[1,1,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,0.01];
ub=[3e3,1e3,1,10,0.1,2,1e2,5,1];

% Define constraints
A=[0 0 0 1 1 -1 0 0 0;
   -1 1 0 0 0 0 0 0 0];
b=[0 0]';

% Initial point
theta0=[415,199,0.16,0.0002,0.0002,0.2,31,4,0.9]';

% Run optimization function
[theta,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta) F3(theta,t_common,Data_median),theta0,A,b,[],[],lb,ub,[],options);

% Compute prediction
t=linspace(t_common(1),t_common(end)*2,1000);
fit = F(theta,t);

% Figures -----------------------------------------------------------------

% Data + prediction
figure('Position', [10 10 800 350])
yyaxis left
plot(t_common,TCELL_median,'*','linewidth',3);
hold on
plot(t,fit(:,1),'-','linewidth',3)
hold off
ylabel('ALC (Cells/uL)')
yyaxis right 
plot(t_common,CART_median,'*','linewidth',3);
hold on
plot(t,fit(:,2),'-','linewidth',3)
ylabel('CAR-T (Cells/uL)')
xlabel('Time (month)')
xlim([0,180])
xticks([0,14,30,60,90,180])
set(gca,'fontsize',20)
%exportgraphics(gcf,'median_fit.png','Resolution',300)


% CART growth rate (rC)
figure()
r=theta(4)+(theta(6)*(fit(:,1)+fit(:,2)-theta(1)).^2)./(theta(5)*(fit(:,1)+fit(:,2)).^2+(fit(:,1)+fit(:,2)-theta(1)).^2);
plot(t,r,'linewidth',3)
set(gca,'fontsize',20)
xticks([0,14,30,60,90,180])
ylabel('Growth rate (day^-1)')
xlabel('Time (month)')
xlim([0,180])
%exportgraphics(gcf,'median_rc.png','Resolution',300)

% CART in log scale
figure()
semilogy(t,fit(:,2),'linewidth',3,'color',[0.8500 0.3250 0.0980])
hold on
semilogy(t_common,Data_median(:,2),'*','linewidth',3);
ylabel('CAR-T (Cells/uL)')
xlabel('Time (month)')
xlim([0,180])
set(gca,'fontsize',20)
%exportgraphics(gcf,'median_logcar.png','Resolution',300)


% CART derivative
figure()
Cdot=-r.*log((fit(:,1)+fit(:,2))/theta(2));
plot(t,Cdot,'linewidth',3)
set(gca,'fontsize',20)
xticks([0,14,30,60,90,180])
ylabel('Proliferation rate (day^-1)')
xlabel('Time (month)')
xlim([0,180])
%exportgraphics(gcf,'median_cdot.png','Resolution',300)


%% Fit individual patients

% Fmincon options
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'Display','none','MaxFunctionEvaluations',9000);

% Bounds and constraints
lb=[1,1,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,0.01];
ub=[3e3,1e3,1,10,0.1,2,1e2,5,1];
A=[0 0 0 1 1 -1 0 0 0;
   -1 1 0 0 0 0 0 0 0];
b=[0 0]';

% Initial estimation
theta0=[269,249,0.2,0.74,1e-4,0.34,30,4.9,0.8]';

% Initialize matrix to store results
Fit_result=zeros(length(IDs),13);

% Loop through all patients

for patient=1:length(IDs)

% Get ALC and CART data
ALC_i=ALC_sub(patient,:);
CART_i=CART_sub(patient,:);

% Remove nans
car_idx=find(~isnan(CART_i));
alc_idx=find(~isnan(ALC_i));
idx=intersect(car_idx,alc_idx);

% Extract Tcell and CART
t_i=t_common(idx);
CART_i=CART_i(idx);
TCELL_i=ALC_i(idx)-CART_i;


% Run Fmincon
[theta,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta) F3(theta,t_i,[TCELL_i',CART_i']),theta0,A,b,[],[],lb,ub,[],options);

% Store results
Fit_result(patient,:)=[theta' fval exitflag, output.firstorderopt, output.stepsize];

% Figures -----------------------------------------------------------------

% Compute prediction
t=linspace(t_common(1),t_common(end),500);
fit = F(theta,t);


% Data + prediction
figure('Position', [10 10 800 350],'visible','off')
yyaxis left
plot(t_i,TCELL_i,'*','linewidth',3);
hold on
plot(t,fit(:,1),'-','linewidth',3)
hold off
ylabel('ALC (Cells/uL)') 
yyaxis right 
plot(t_i,CART_i,'*','linewidth',3);
hold on
plot(t,fit(:,2),'-','linewidth',3)
ylabel('CAR-T (Cells/uL)')
xlabel('Time (month)')
xlim([0,180])
ylim([0,150])
title(['patient ',IDs{patient}])
xticks([0,14,30,60,90,180])
set(gca,'fontsize',20)
%exportgraphics(gcf,[IDs{patient} '_fit.png'],'Resolution',300)


% CART growth rate (rC)
figure('visible','off')
r=theta(4)+(theta(6)*(fit(:,1)+fit(:,2)-theta(1)).^2)./(theta(5)*(fit(:,1)+fit(:,2)).^2+(fit(:,1)+fit(:,2)-theta(1)).^2);
plot(t,r,'linewidth',3)
set(gca,'fontsize',20)
xticks([0,14,30,60,90,180])
ylabel('Proliferation rate (day^-1)')
xlabel('Time (month)')
xlim([0,180])
%exportgraphics(gcf,[IDs{patient} '_rc.png'],'Resolution',300)


% CART in log scale
figure('visible','off')
semilogy(t,fit(:,2),'linewidth',3,'color',[0.8500 0.3250 0.0980])
hold on
semilogy(t_i,CART_i,'*','linewidth',3);
ylabel('CAR-T (Cells/uL)')
xlabel('Time (month)')
xlim([0,180])
set(gca,'fontsize',20)
%exportgraphics(gcf,[IDs{patient} '_logcar.png'],'Resolution',300)


% CART derivative
figure('visible','off')
Cdot=-r.*log((fit(:,1)+fit(:,2))/theta(2));
plot(t,Cdot,'linewidth',3)
set(gca,'fontsize',20)
xticks([0,14,30,60,90,180])
ylabel('Proliferation rate (day^-1)')
xlabel('Time (month)')
xlim([0,180])
%exportgraphics(gcf,[IDs{patient} '_cdot.png'],'Resolution',300)

end



%% Fit every patient with initial estimation sampling

% Define bounds and constraints
lb=[1,1,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,0.01];
ub=[3e3,1e3,1,10,0.1,2,1e2,5,1];
A=[0 0 0 1 1 -1 0 0 0];
b=0;

% Define number of sampled initial estimations
N=2;

% initialize matrix to store results
Fit_result2=zeros(length(IDs),13);


% Parallel loop
parfor patient=1:length(IDs)

% Get Tcell and CART values

ALC_i=ALC_sub(patient,:);
CART_i=CART_sub(patient,:);

car_idx=find(~isnan(CART_i));
alc_idx=find(~isnan(ALC_i));
idx=intersect(car_idx,alc_idx);

t_i=t_common(idx);
CART_i=CART_i(idx);
TCELL_i=ALC_i(idx)-CART_i;

% Quasi-random sample 9 dimensional space (Sobol)

p=sobolset(9);
p=scramble(p,'MatousekAffineOwen');
param=net(p,N);

% Transform [0,1] to intervals for each parameter [10^a,10^b]

kN_int=10.^(param(:,1)+2);
kC_int=10.^(param(:,2)+2);
rN_int=10.^(param(:,3)*2-2);
pC_int=10.^(param(:,4)*4-4);
a_int=10.^(param(:,5)*4-5);
b_int=10.^(param(:,6)*5-5);
N0_int=10.^(param(:,7)*5-3);
C0_int=10.^(param(:,8)*3-3);
c_int=10.^(param(:,9)*2-2);

% Initialize vector to store residual errors
error=[];

% Loop through all initial estimations
for i=1:N

% Extract initial estimation
theta0=[kN_int(i),kC_int(i),rN_int(i),pC_int(i),a_int(i),b_int(i),N0_int(i),C0_int(i),c_int(i)]';

% Ensure parameter bounds are not exceeded
for j=1:9
    if theta0(j)<lb(j)
        theta0(j)=1.000001*lb(j);
    elseif theta0(j)>ub(j)
        theta0(j)=0.999999*ub(j);
    end
end

% Run optimization 
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',800);
[theta,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta) F3(theta,t_i,[TCELL_i',CART_i']),theta0,A,b,[],[],lb,ub,[],options);

% Store error
error(i)=fval;
end

% Select estimation with minimum error and store 
[val,idx]=min(error);
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'Display','iter','MaxFunctionEvaluations',9000);
theta0=[kN_int(idx),kC_int(idx),rN_int(idx),pC_int(idx),a_int(idx),b_int(idx),N0_int(idx),C0_int(idx),c_int(idx)]';
[theta,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta) F3(theta,t_i,[TCELL_i',CART_i']),theta0,A,b,[],[],lb,ub,[],options);
Fit_result2(patient,:)=[theta' fval exitflag, output.firstorderopt,output.stepsize];

end

% Save estimation

writematrix(Fit_result2,"Parameters_NC.xlsx")



%% Profile likelihood computation 

% Define options and constraints
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'Display','iter','MaxFunctionEvaluations',9000);
A=[0 0 0 1 1 -1 0 0 0;
   -1 1 0 0 0 0 0 0 0];
b=[0 0]';

% Number of points in interval for each parameter
n=10;

% Number of initial estimations
nruns=200;

% Standard deviation for initial estimation sampling
sigma=0.05;

% Initialize structures to store results
INT=zeros(length(IDs),9,n);
PLIK=zeros(length(IDs),9,n);

% Loop through patients

parfor patient=1:length(IDs)

% Extract data

ALC_i=ALC_sub(patient,:);    
CART_i=CART_sub(patient,:);

car_idx=find(~isnan(CART_i));
alc_idx=find(~isnan(ALC_i));
idx=intersect(car_idx,alc_idx);

t_i=t_common(idx);
CART_i=CART_i(idx);
TCELL_i=ALC_i(idx)-CART_i; 

% Define bounds

lb=[1,1,1e-5,1e-5,1e-5,1e-5,1e-3,1e-3,0.01];
ub=[3e3,1e3,1,10,0.1,2,1e2,5,1];

% Extract parameters from best fit
fit=Fit_result(patient,1:9);

% Loop through all parameters

for param=1:9

% Define interval for parameter
int=linspace(0.1*fit(param),1.9*fit(param),n);

% Initialize variables to store results (likelihood and residuals)

plik=zeros(nruns,n);
fvals=zeros(1,nruns);

% Loop through points in interval

for i=1:length(int)

    % Do nruns estimations and select best

    for k=1:nruns
    
    % Sample initial estimation except for selected parameter
    theta0=normrnd(fit,sigma*fit);
    theta0(param)=int(i);
    lb(param)=int(i);
    ub(param)=int(i);

    % Ensure bounds are not exceeded

    for j=1:9
        if theta0(j)<lb(j)
            theta0(j)=1.000001*lb(j);
        elseif theta0(j)>ub(j)
            theta0(j)=0.999999*ub(j);
        end
    end
 
    % Run optimization function

    [theta,fval] = fmincon(@(theta) F3(theta,t_i,[TCELL_i',CART_i']),theta0,A,b,[],[],lb,ub,[],options);

    % Store residual
    plik(k,i)=fval;

    end

end

% Store results

INT(patient,param,:)=int;
PLIK(patient,param,:)=min(plik);

end

end


%% Functions

% Compute predictions for a given parameter estimation

function C=F(theta,t)

c0=[theta(7);theta(8)];
%c0=[6,0.36];
[time,C]=ode45(@model,t,c0);

    function dndt=model(t,n)

    dndt=zeros(2,1);

    kN=theta(1);    
    kC=theta(2);   
    rN=theta(3);      
    pC=theta(4);   
    a=theta(5);    
    b=theta(6);
    c=theta(9);
 
    T=n(1)+n(2);
    rc=pC+b*((T-kN)^2)/(a*T^2+(T-kN)^2);
    dndt(1)=-rN*n(1)*log((n(1)+c*n(2))/kN);
    dndt(2)=-rc*n(2)*log(T/kC);

    end

end

% Optimization function

function ssr=F3(theta,xdata,ydata)

c0=[theta(7);theta(8)];

[time,C]=ode45(@model,xdata,c0);

alcdata=ydata(1:end,1);
cardata=ydata(2:end,2);
alcest=C(1:end,1);
carest=C(2:end,2);

ssr_alc=sum((log(abs(alcest))-log(abs(alcdata))).^2);
ssr_car=sum((log(abs(carest))-log(abs(cardata))).^2);
ssr=1*ssr_alc+1*ssr_car;


    function dndt=model(t,n)

    dndt=zeros(2,1);

    kN=theta(1);    
    kC=theta(2);   
    rN=theta(3);      
    pC=theta(4);   
    a=theta(5);    
    b=theta(6);
    c=theta(9);

    T=n(1)+n(2);
    rc=pC+b*((T-kN)^2)/(a*T^2+(T-kN)^2);
    dndt(1)=-rN*n(1)*log((n(1)+c*n(2))/kN);
    dndt(2)=-rc*n(2)*log(T/kC);

    if n(2)<1e-6
        dndt(2)=0;
    end

    end
end





