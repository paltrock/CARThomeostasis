%% Data

Data=readtable("/Users/alvaro/Desktop/Full_dataset.xlsx");
Data([4,15,17,18],:)=[]; % Remove KTE 6,14,20,21 
ALC=Data{:,18:25}*1000*5e6; % 5L=5e6uL (ref?)
CART=Data{:,26:36}*5e6; % 5L=5e6uL (ref?)
Tumor_initial=Data{:,14}*1e8; % 1mL=1cm3=1e8 cells (del monte 2009, wikipedia)
Tumor_response=readmatrix("/Users/alvaro/Desktop/Tumor_response.xlsx");
Tumor_response=Tumor_response(:,2:4);
IDs=Data{:,1}; % Patients' IDs

Parameters=readmatrix("/Users/alvaro/Desktop/Code/Parameters_NC.xlsx");
Parameters=readmatrix("/Users/alvaro/Desktop/Full_parameters.xlsx");
Parameters=Parameters(:,2:end); 
%Parameters([4,15,17,18],:)=[]; % Remove KTE 6,14,20,21 
Parameters(:,[1:2,7:8])=Parameters(:,[1:2,7:8])*5e6; 

tint=0.1;
tspan=0:tint:180;


%% Range for rB


% Size of intervals
grid=15; % kB and gB
rB_size=50; % rB

% Specify intervals for the three parameters
rB=logspace(-3,0,rB_size);
gB=logspace(0,2,grid);
kB=logspace(9,11,grid);

% Initialize variables to store results
rB_min=zeros(length(IDs),length(rB));
rB_range=zeros(length(IDs),3);

% Loop all patients
parfor patient=1:length(IDs)

% Initialize residual matrix
res=zeros(grid);

% For each rB
for k=1:rB_size

    % For each pair (kB,gB)
for i=1:grid
    for j=1:grid
        
        % Select parameters
        param=[Parameters(patient,[1:6,9]) rB(k) gB(i) kB(j)];
        
        % Set initial state
        N0=Parameters(patient,7);
        C0=Parameters(patient,8);
        B0=Tumor_initial(patient);
        n0=[N0,C0,B0];

        % Solve model
        Opt    = odeset('Events', @myEvent);
        [t,n] = ode45(@(t,n) model(t,n,param),tspan,n0',Opt);

        % Extract tumor prediction
        B=n(:,3);

        % Compute response at days 30, 90, 180 and measure residual

        R=response(B,t);

        res(i,j)=sum(abs(R-Tumor_response(patient,:)));

    end
end

% Find minimum rB in the whole grid

rB_min(patient,k)=min(min(res));

end

end

%% Figure

for patient=1:length(IDs)

figure()
semilogx(rB,rB_min(patient,:),'LineWidth',3)
xlabel('rB')
ylabel('minimum deviation')
set(gca,'fontsize',15)
%exportgraphics(gcf,[IDs{patient},'_rB.png'],'Resolution',300)

% Compute range of minimum for rB and store

idx=find(rB_min(patient,:)==min(rB_min(patient,:)));
rB_range(patient,:)=[min(rB(idx)),max(rB(idx)),min(rB_min(patient,:))];

end

rB_data=max(rB_range');

%% Fit gB and kB (fixed rB)

% Define intervals
grid=50;
gB=logspace(0,2,grid);
kB=logspace(9,11,grid);


% Initialize variables
gB_range=zeros(length(IDs),2);
kB_range=zeros(length(IDs),2);
gB_selectionA=zeros(length(IDs),1);
kB_selectionA=zeros(length(IDs),1);
gB_selectionB=zeros(length(IDs),1);
kB_selectionB=zeros(length(IDs),1);

% Loop through all patients
parfor patient=1:length(IDs)

res=zeros(grid);

for i=1:grid
    for j=1:grid

        param=[Parameters(patient,[1:6,9]) rB_data(patient) gB(i) kB(j)];

        N0=Parameters(patient,7);
        C0=Parameters(patient,8);
        B0=Tumor_initial(patient);
        n0=[N0,C0,B0];

        Opt    = odeset('Events', @myEvent);
        [t,n] = ode45(@(t,n) model(t,n,param),tspan,n0',Opt);

        B=n(:,3);

        % Compute response at days 30, 90, 180 and measure error

        R=response(B,t);

        res(i,j)=sum(abs(R-Tumor_response(patient,:)));

    end
end

% Obtain optimum gB,kB (see notes)

% Get all pairs for which res is minimum
[gB_idx, kB_idx]=find(res==min(min(res)));

% Get ranges
kB_range(patient,:)=[min(kB(kB_idx)),max(kB(kB_idx))];
gB_range(patient,:)=[min(gB(gB_idx)),max(gB(gB_idx))];

% Option A (median kB, corresponding gB)

median_kBA=kB_idx(floor(median(1:length(unique(kB_idx)))));
median_gBA=floor(median(gB_idx(kB_idx==median_kBA)));

kB_selectionA(patient)=kB(median_kBA);
gB_selectionA(patient)=gB(median_gBA);

% Option B (median gB, corresponding kB)

median_gBB=gB_idx(floor(median(1:length(unique(gB_idx)))));
median_kBB=floor(median(kB_idx(gB_idx==median_gBB)));

kB_selectionB(patient)=kB(median_kBB);
gB_selectionB(patient)=gB(median_gBB);

% Figure ------------------------------------------------------------

% Heatmap

figure()
imagesc(kB,gB,res)
set(gca,'xscale','log','yscale','log','Ydir','normal')
set(gca,'fontsize',10,'linewidth',2)
colormap(parula(10));
colorbar;
xlabel('kB')
ylabel('gB')
clim([0,10])
xlim([min(kB),max(kB)])
ylim([min(gB),max(gB)])
hold on

% Selected points

plot(kB(median_kBA),gB(median_gBA),'r+','markersize',20,'linewidth',3)
plot(kB(median_kBB),gB(median_gBB),'y+','markersize',20,'linewidth',3)

text(kB(median_kBB),gB(median_gBB),'B','fontsize',20,'color','yellow')
text(kB(median_kBA),gB(median_gBA),'A','fontsize',20,'color','red')


%exportgraphics(gcf,[IDs{patient},'_hmselection.png'],'Resolution',300)


end

gB_data=median([gB_selectionA,gB_selectionB],2);
kB_data=median([kB_selectionA,kB_selectionB],2);

% Save tumor parameters

writematrix([rB_data',gB_data,kB_data],"tumonr_parameters.xlsx")




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

% Function to stop solver when B>1e12
function [value, isterminal, direction] = myEvent(T, Y)
value      = (Y(3) > 1e12);
isterminal = 1;   % Stop the integration
direction  = 0;
end

% Compute response at days 30,90,180 from tumor data

function R=response(B,t)

tint=t(2)-t(1);
R=zeros(1,3);
T=length(t);

        idx30=30/tint;
        if idx30<T
 
            if B(idx30)<1e8
                R(1)=0;
            elseif B(idx30)<B(1)
                R(1)=1;
            elseif B(idx30)>B(1)
                R(1)=2;
            end

        else
            R(1)=3;
        end

    

        idx90=90/tint;
        if idx90<T

            if B(idx90)<1e8
                R(2)=0;
            elseif B(idx90)<B(idx30)
                R(2)=1;
            elseif B(idx90)>B(idx30)
                R(2)=2;
            end

        else
            R(2)=3;
        end



        idx180=180/tint;
        if idx180<T

            if B(idx180)<1e8
                R(3)=0;
            elseif B(idx180)<B(idx90)
                R(3)=1;
            elseif B(idx180)>B(idx90)
                R(3)=2;
            end

        else
            R(3)=3;
        end
        


end