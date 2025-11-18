% ESTIMATE HYBRID MODEL FOR VARIOUS EPSILON

 %% Set paths

my_dir = "/Users/alvaromartinezrubio/Desktop/individualCAR-T-main"; 
cd(my_dir)
addpath(genpath(my_dir))

%% Import data

% Read dataset
Data=readtable("Data_full.xlsx");

% Retrieve T-cell, CAR T-cell and tumor data
CART=Data{:,2:12};
TCELL=Data{:,13:23};
TUM=Data{:,24:34};

% Define measurement precision (for parameter estimation)
TCELL_sd=10*ones(1,size(TCELL,2));
CART_sd=2e-3*ones(1,size(CART,2));
TUM_sd=0.2*ones(1,size(TUM,2));

% Retrieve patient's IDs
IDs=Data{:,1};

% Time points
time=[0 7 14 30 90 180 270 360 450 540 720];

% Parameter estimates
Fit_H=readmatrix("Parameters_H.xlsx");
Fit_A=readmatrix("Parameters_A.xlsx");


%% Estimate hybrid model

% Define number of sampled initial estimations

Np=100;

% initialize matrix to store parameter estimates

Fit=zeros(length(IDs),18);


% Construct LB and UB matrices
%  [kN ,kC , rN , pC , aH ,aA , b  , rB , gB ,kB , eps, N0 , C0 , B0 ]
%  [kN ,kC , rN , pC , a , b  , rB , gB ,kB , N0 , C0 , B0  ]

LB=zeros(length(IDs),14);
UB=zeros(length(IDs),14);
iHA=[1,2,3,4,7,8,9,10,12,13,14];
iHyA=[1,2,3,4,6,7,8,9,10,11,12];

for i=1:length(iHA)

LB(:,iHA(i))=0.5*min(Fit_A(:,iHyA(i)),Fit_H(:,iHyA(i)));
UB(:,iHA(i))=1.5*max(Fit_A(:,iHyA(i)),Fit_H(:,iHyA(i)));

end

LB(:,5)=Fit_H(:,5);
UB(:,5)=Fit_H(:,5);
LB(:,6)=Fit_A(:,5);
UB(:,6)=Fit_A(:,5);

% Convert to cell for parfor
LBcell=num2cell(LB,2);
UBcell=num2cell(UB,2);

% Constraints

A=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0]; 
b=0;

% Epsilon

eps_values=[0,0.25,0.5,0.75,1];
eps_SSR=zeros(length(IDs),length(eps_values));
eps_SSRcell=num2cell(eps_SSR,2);
Fit_aux = cell(1, length(IDs));
Fit_aux(:) = {zeros(length(eps_values),18)};


% Parallel loop

baseFile = ('part_%03d.mat');
tic
parfor patient=1:length(IDs)

% Get ALC and CART data
TCELL_i=TCELL(patient,:);
CART_i=CART(patient,:);
TUM_i=TUM(patient,:);

% Remove nans
car_idx=find(~isnan(CART_i));
tcell_idx=find(~isnan(TCELL_i));
tum_idx=find(~isnan(TUM_i));

% Extract Tcell and CART
CART_t_i=time(car_idx);
TCELL_t_i=time(tcell_idx);
TUM_t_i=time(tum_idx);
CART_i=CART_i(car_idx);
TCELL_i=TCELL_i(tcell_idx);
TUM_i=TUM_i(tum_idx);
CART_sd_i=CART_sd(car_idx);
TCELL_sd_i=TCELL_sd(tcell_idx);
TUM_sd_i=TUM_sd(tum_idx);

N=[TCELL_t_i' TCELL_i' TCELL_sd_i'];
C=[CART_t_i' CART_i' CART_sd_i'];
B=[TUM_t_i' TUM_i' TUM_sd_i'];

for j=1:length(eps_values)

% Fix epsilon

LBcell{patient}(11)=eps_values(j);
UBcell{patient}(11)=eps_values(j);

% Quasi-random sample parameter space (Sobol)

p=sobolset(length(A));
p=scramble(p,'MatousekAffineOwen');
p=net(p,Np);

% Transform [0,1] to intervals for each parameter

param=10.^(p.*(log10(UBcell{patient})-log10(LBcell{patient}))+log10(LBcell{patient}));

% Initialize vector to store residual errors
error=ones(1,Np);

% Loop through all initial estimations
for i=1:Np

% Extract initial estimation
theta0=param(i,:)';

% Run optimization 
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'MaxFunctionEvaluations',9000,'display','off');
[~,fval] = fmincon(@(theta) SSR_HA(theta,N,C,B),theta0',A,b,[],[],LBcell{patient},UBcell{patient},[],options);

% Store error
error(i)=fval;
end

% Select estimation with minimum error and store 
[val,idx]=min(error);
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'Display','off','MaxFunctionEvaluations',9000);
theta0=param(idx,:)';
[theta,fval,exitflag,output] = fmincon(@(theta) SSR_HA(theta,N,C,B),theta0',A,b,[],[],LBcell{patient},UBcell{patient},[],options);
Fit_aux{patient}(j,:)=[theta fval exitflag, output.firstorderopt,output.stepsize];

eps_SSRcell{patient}(j)=fval;

end

[val,idx]=min(eps_SSRcell{patient});
singleRecord=Fit_aux{patient}(idx,:);
Fit(patient,:)=singleRecord;

parsave(sprintf(baseFile, patient), singleRecord);

end

eps_SSR=cell2mat(eps_SSRcell);
toc
writematrix(Fit,"Parameters_eps.xlsx")

%% Load and combine all temporary files

Fit_temp=zeros(length(IDs),18);
baseFile = ('part_%03d.mat');
for i = 1:length(IDs)
    try
    d = load(sprintf(baseFile, i));
    Fit_temp(i,:) = d.var;
    catch ME
        ME
    end
end

writematrix(Fit_temp,"Parameters_eps_Temp.xlsx")

%% Functions

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Compute predictions for a given parameter estimation

function C=F_HA(theta,t)

c0=[theta(12);theta(13);theta(14)];


options = odeset('Events',@myEvent);
[time,C]=ode45(@model,t,c0,options);


    function dndt=model(t,n)

    dndt=zeros(3,1);

    kN=theta(1);    
    kC=theta(2);   
    rN=theta(3);      
    pC=theta(4);   
    aH=theta(5);    
    aA=theta(6);
    b=theta(7);
    rB=theta(8);
    gB=theta(9);
    kB=theta(10);
    eps=theta(11);

    T=n(1)+n(2);
    rc=pC+b*(eps*((T-kN)^2)/(aH*T^2+(T-kN)^2)+(1-eps)*(n(3)^2)/(n(3)^2+aA));
    dndt(1)=-rN*n(1)*log(T/kN);
    dndt(2)=-rc*n(2)*log(T/kC);
    dndt(3)=rB*n(3)-gB*n(3)*n(2)/(kB+n(2));

    if n(1)<2e-7
        dndt(1)=0;
    end

    if n(2)<2e-7
        dndt(2)=0;
    end

    if n(3)<2e-7
        dndt(3)=0;
    end

    if n(3)>1e12 % Stop if tumor grows to max size
        dndt(1)=0;
        dndt(2)=0;
        dndt(3)=0;
    end

    end

end

% -------------------------------------------------------------------------


% Optimization function

function ssr=SSR_HA(theta,N,C,B)

c0=[theta(12);theta(13);theta(14)];

t_end=B(end,1);

t_int=0.1;

options = odeset('Events',@myEvent);
[time,y]=ode45(@model,0:t_int:t_end,c0,options);

if time(end)<t_end

    ssr=1e20;

else

Ndata=N(:,2);
Cdata=C(:,2);
Bdata=B(:,2);
Nest=y(N(:,1)/t_int+1,1);
Cest=y(C(:,1)/t_int+1,2);
Best=y(B(:,1)/t_int+1,3);
Nsd=N(:,3);
Csd=C(:,3);
Bsd=B(:,3);

ssr_alc=sum(((Nest-Ndata)./(Nsd.*Ndata)).^2);
ssr_car=sum(((Cest-Cdata)./(Csd.*Cdata)).^2);
ssr_tum=sum(((Best-Bdata)./(Bsd.*Bdata)).^2);

ssr=1*ssr_alc+0.000001*ssr_car+0.001*ssr_tum; 

end


    function dndt=model(t,n)

    dndt=zeros(3,1);

    kN=theta(1);    
    kC=theta(2);   
    rN=theta(3);      
    pC=theta(4);   
    aH=theta(5);    
    aA=theta(6);
    b=theta(7);
    rB=theta(8);
    gB=theta(9);
    kB=theta(10);
    eps=theta(11);

    T=n(1)+n(2);
    rc=pC+b*(eps*((T-kN)^2)/(aH*T^2+(T-kN)^2)+(1-eps)*(n(3)^2)/(n(3)^2+aA));
    dndt(1)=-rN*n(1)*log(T/kN);
    dndt(2)=-rc*n(2)*log(T/kC);
    dndt(3)=rB*n(3)-gB*n(3)*n(2)/(kB+n(2));

    if n(1)<2e-7
        dndt(1)=0;
    end


    if n(2)<2e-7
        dndt(2)=0;
    end

    if n(3)<2e-7
        dndt(3)=0;
    end

    if n(3)>1e12 % Stop if tumor grows to max size
        dndt(1)=0;
        dndt(2)=0;
        dndt(3)=0;
    end

    end
end

% -------------------------------------------------------------------------


% Function to stop solver when B>1e13 cells
function [value, isterminal, direction] = myEvent(T, Y)
value      = (Y(3) > 1e13/5e6);
isterminal = 1;   % Stop the integration
direction  = 0;
end

% -------------------------------------------------------------------------

function parsave(fname, var)
  save(fname, 'var')
end
