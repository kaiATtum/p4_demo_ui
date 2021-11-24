clear all;
clc;

%% get psp price 
psp_price = readtable('MEP_02-Sep-2018_to_02-Sep-2018.csv');
load_profilex = xlsread('20180827.xls');

number_of_market_period = size(psp_price,1);

%number_of_market_period = number_of_market_period/6;

number_of_market_period = 12;

psp_pricer = table2array(psp_price(:,'PRICE___MWh_'));

psp_price = zeros(number_of_market_period,1);
load_profile = load_profilex(1:2:2*number_of_market_period,19)/max(load_profilex(:,19)); % data on 02/09

for i = 1:number_of_market_period
    psp_price(i) = psp_pricer(2*i-1);
end


%%
%mpc = case33_radial_test;
mpc = case141_with_storageV3;
define_constants;

%% Verification of linDist FLow 
%Ybus = makeYbus(mpc);
[Ybus1, Yf, Yt] = makeYbus(mpc);
f_bus_idx = mpc.branch(:,F_BUS);
t_bus_idx = mpc.branch(:,T_BUS);

nbr = size(mpc.branch,1);
nbu = size(mpc.bus,1);
nstroage = size(mpc.storage,1);

L = zeros(nbu,nbu);
ngen = size(mpc.gen,1);

for br = 1:nbr
    br_F_BUS = mpc.branch(br,F_BUS);
    br_T_BUS = mpc.branch(br,T_BUS);
    br_BR_R = mpc.branch(br,BR_R);
    br_BR_X = mpc.branch(br,BR_X);
    br_Y = 1 / (br_BR_R + 1j * br_BR_X);
    
    L(br_F_BUS, br_T_BUS) = -br_Y;
    L(br_T_BUS, br_F_BUS) = -br_Y;
    L(br_F_BUS, br_F_BUS) = L(br_F_BUS, br_F_BUS) + br_Y;
    L(br_T_BUS, br_T_BUS) = L(br_T_BUS, br_T_BUS) + br_Y;
end
Ybus = L; 

ReY = real(Ybus(2:end,2:end));
ImY = imag(Ybus(2:end,2:end));

%Rbus = inv(ReY);
%Xbus = inv(ImY);

%v_lin = ones(nbu,1)+ Rbus*mpc.bus(:,PD)/mpc.baseMVA + Xbus*mpc.bus(:,QD)/mpc.baseMVA;

LinDist = inv([ReY -ImY; -ImY -ReY]); 
vbase = mpc.gen(1,VG)*ones(nbu-1,1);
u_lin = LinDist*[mpc.bus(2:end,PD)/mpc.baseMVA; mpc.bus(2:end,QD)/mpc.baseMVA];
v_lin = -u_lin(1:(nbu-1)) + vbase; 
a_lin = -u_lin(nbu:end);

% voltage - good 
pfres = runpf(mpc);
v_matpower =  pfres.bus(:,VM);

% figure()
% plot(v_matpower(2:end),'r-*'); hold on;
% plot(v_lin,'b-o'); hold on;

% angle - bad - angles are not in the power flow formulation 
% a_matpower =  pfres.bus(:,VA)/180*2*pi;
% figure()
% plot(a_matpower(2:end),'r-*'); hold on;
% plot(a_lin,'b-o'); hold off;

% loss - bad - losless estimation 
loss_linDist = mpc.gen(1,VG)*ones(1,nbu-1)*Ybus(2:end,2:end)*mpc.gen(1,VG)*ones(nbu-1,1);
loss = get_losses(pfres);


% Z01 = 1/Ybus(1,1);
% Z12 = (mpc.branch(2,BR_R) + mpc.branch(2,BR_X)*1i)/(Vbase^2 / Sbase);
% V0 = mpc.gen(1,VG);
% V1 = v_lin(1);
% V2 = v_lin(2);
% S_Load1 = mpc.bus(2,PD) +  mpc.bus(2,QD)*1i;

S0 = sum(mpc.bus(:,PD))/mpc.baseMVA + sum(mpc.bus(:,QD))*1i/mpc.baseMVA;
P0 = sum(mpc.bus(:,PD))/mpc.baseMVA;
Q0 = sum(mpc.bus(:,QD))/mpc.baseMVA;

%% %% LinDistOPF centralized solution 
time_step = number_of_market_period; %market period 0.5h 
soc_init = [0.3 0.6 0.55]';
Bess_capacity = mpc.storage(:,2);
Bess_charge_efficiency = mpc.storage(:,3);
Bess_discharge_efficiency = mpc.storage(:,4);
Bess_cost_coefficient = mpc.storage(:,6);


QG_max = zeros(size(mpc.bus,1),1);
PG_max = zeros(size(mpc.bus,1),1);
QG_min = zeros(size(mpc.bus,1),1);
PG_min = zeros(size(mpc.bus,1),1);

QG_max(mpc.gen(:,1)) = mpc.gen(:,QMAX)/mpc.baseMVA;
PG_max(mpc.gen(:,1)) = mpc.gen(:,PMAX)/mpc.baseMVA;
QG_min(mpc.gen(:,1)) = mpc.gen(:,QMIN)/mpc.baseMVA;
PG_min(mpc.gen(:,1)) = mpc.gen(:,PMIN)/mpc.baseMVA;

% ramping rate limit of energy storage 
PG_max(mpc.storage(:,1)) = mpc.storage(:,5)/mpc.baseMVA;
PG_min(mpc.storage(:,1)) = -mpc.storage(:,5)/mpc.baseMVA;


% Generation Set points
PG_var = sdpvar(size(mpc.bus,1),time_step);
QG_var = sdpvar(size(mpc.bus,1),time_step);

P_charge_var = sdpvar(size(mpc.storage,1),time_step);
P_discharge_var = sdpvar(size(mpc.storage,1),time_step);
SOC_var = sdpvar(size(mpc.storage,1),time_step);

PQ_inj_var = [PG_var(2:end,:) - repmat(mpc.bus(2:end,PD)/mpc.baseMVA,1,time_step)*diag(load_profile); QG_var(2:end,:) - repmat(mpc.bus(2:end,QD)/mpc.baseMVA,1,time_step)*diag(load_profile)];
deltaVA = - LinDist*PQ_inj_var; % assume the load does not change 

%% robust optimization varibales 
% robust price
    mu_sum = sdpvar(1,1);
    mu_prox = sdpvar(1,time_step);

    % parameters for robust optimization 
    level_of_conserv = 3;
    max_deviation_from_psp_price = 0.1*psp_price;

% robust load
    level_of_conserv_D = 0.04;  
    mu_sum_D = sdpvar(1,1);
    mu_prox_D = sdpvar(1,time_step);

    max_deviation_load = 0.1*load_profile;

%% constraint time step loop
PQ0_constr = [];
V_constr = [];
PQ_constr = [];
Batt_constr = [];
SOC_constr = [];
robust_constr = [];

SOC_min = 0.3;
SOC_max = 1;

for t = 1:time_step
% p_0, q_0 constr
    PQ0_constr = [PQ0_constr, PG_var(1,t) == P0*load_profile(t)-sum(PG_var(2:end,t)) + mu_prox_D(t) + level_of_conserv_D*mu_sum_D, QG_var(1,t) == Q0*load_profile(t)-sum(QG_var(2:end,t))];
    V_constr = [V_constr, 0.90*ones(nbu-1,1)/1.001 <= deltaVA(1:(nbu-1),t) + vbase <= 1.05*ones(nbu-1,1)*1.00]; % 5% maximal error is tolerated
    PQ_constr = [PQ_constr, PG_min(2:end)<=PG_var(2:end,t)<=PG_max(2:end),QG_min(2:end)<=QG_var(2:end,t)<=QG_max(2:end), PG_min(1)<= PG_var(1,t) <= PG_max(1), QG_min(1)<= QG_var(1,t)  <= QG_max(1)];
    Batt_constr = [Batt_constr, PG_var(mpc.storage(:,1),t) == P_discharge_var(:,t) - P_charge_var(:,t), P_discharge_var(:,t)>=0, P_charge_var(:,t) >= 0];
    
    if t == 1
        SOC_constr = [SOC_var(:,t) == soc_init - inv(diag(Bess_capacity))*(inv(diag(Bess_discharge_efficiency))*P_discharge_var(:,t) - inv(diag(Bess_charge_efficiency))*P_charge_var(:,t)), SOC_var(:,t) >= SOC_min, SOC_var(:,t) <= SOC_max];
    else
        SOC_constr = [SOC_constr, SOC_var(:,t) == SOC_var(:,t-1) - inv(diag(Bess_capacity))*(inv(diag(Bess_discharge_efficiency))*P_discharge_var(:,t) - inv(diag(Bess_charge_efficiency))*P_charge_var(:,t)), SOC_var(:,t)>=SOC_min, SOC_var(:,t)<=SOC_max];
    end
end

% robust constraint for price 
robust_constr = [mu_sum>=0];

for t = 1:time_step
    robust_constr = [robust_constr, mu_sum + mu_prox(t) >= max_deviation_from_psp_price(t)*PG_var(1,t), mu_prox(t)>= 0]; 
end

% robust constraint for the load
robust_constr = [robust_constr, mu_sum_D>=0];

for t = 1:time_step
    robust_constr = [robust_constr, mu_sum_D + mu_prox_D(t) >= max_deviation_load(t)*P0, mu_prox_D(t)>= 0]; 
end
%% objective 
%assume the price doesnot vary for DG

dg_cost = 0;
battery_cost = 0;
psp_cost = 0;
robust_cost = 0;

for t = 1:time_step
    dg_cost = dg_cost + (mpc.gencost(2:ngen,5)'*PG_var(mpc.gen(2:end,1),t) + PG_var(mpc.gen(2:end,1),t)'*diag(mpc.gencost(2:ngen,6))'*PG_var(mpc.gen(2:end,1),t)  + mpc.gencost(ngen+2:end,5)'*QG_var(mpc.gen(2:end,1),t));
    for c = 1:3
        battery_cost = battery_cost + Bess_cost_coefficient(c)*(P_discharge_var(c,t)-P_charge_var(c,t))^2;
    end
    psp_cost = psp_cost + PG_var(1,t)*psp_price(t) + QG_var(1,t)*psp_price(t)/3;
end

% robustification objective 
for t = 1:time_step
    robust_cost = robust_cost + (mu_prox(t) + PG_var(1,t)*psp_price(t));
end

robust_cost = robust_cost + level_of_conserv*mu_sum;

price_factor = 1;

cost = (price_factor*battery_cost + psp_cost + price_factor*dg_cost) + robust_cost; 
constr = [PQ0_constr,V_constr,PQ_constr,Batt_constr,SOC_constr,robust_constr]; 

ops = sdpsettings('solver','gurobi');

optimize(constr,cost,ops);

soc_res = value(SOC_var);

Pvar_res = value(PG_var);

cost_total_central = value(cost);

cost_robust_dual = value(robust_cost);

v_mag_central_res = value(deltaVA(1:(nbu-1),2) + vbase);

energy_import_psp_res = Pvar_res(1,:)';


%% derive the robust prices - worst case scenario 
PCC_schedule = Pvar_res(1,:)';

u_var = sdpvar(1,time_step);
constr_primal = 0;
obj_primal = 0;

for t = 1:time_step
    obj_primal = obj_primal + PCC_schedule(t)*(psp_price(t) + max_deviation_from_psp_price(t)*u_var(t));
end

for t = 1:time_step
    constr_primal = [constr_primal, u_var(t)<=1, u_var(t) >= 0];
end
constr_primal = [constr_primal, sum(u_var) <= level_of_conserv];
obj_primal = - obj_primal; % maximization problem 


optimize(constr_primal,obj_primal,ops);

cost_robust_primal = value(obj_primal);
u_value = value(u_var);

robust_price = psp_price + max_deviation_from_psp_price.*u_value';

%% %% derive the robust load - worst case scenario 
Robust_load = P0*load_profile + value(mu_prox_D)'+value(mu_sum_D)*level_of_conserv_D;

figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;

plot(Robust_load,'r--*'); hold on;
plot(P0*load_profile,'ko'); hold off;


%% plot results 
figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;

plot(robust_price,'r--*'); hold on;
plot(psp_price+max_deviation_from_psp_price); hold on;
plot(psp_price-max_deviation_from_psp_price); hold on;
plot(psp_price); hold off;
% Qvar_res = value(QG_var);
% 
% PGres = Pvar_res(mpc.gen(:,1));
% QGres = Pvar_res(mpc.gen(:,1));
% 
% %% new linDist Power Flow
% u_lin_new = LinDist*[-Pvar_res(2:end) + mpc.bus(2:end,PD)/mpc.baseMVA; -Qvar_res(2:end) + mpc.bus(2:end,QD)/mpc.baseMVA];
% v_lin_new = -u_lin_new(1:(nbu-1)) + vbase; 
% plot(v_lin_new,'yo'); hold on; 
% 
% mpc.gen(2:end,PG) = Pvar_res(mpc.gen(2:end,1))*mpc.baseMVA;
% mpc.gen(2:end,QG) = Qvar_res(mpc.gen(2:end,1))*mpc.baseMVA;
% 
% pfres_new = runpf(mpc);
% v_matpower_new =  pfres_new.bus(:,VM);
% plot(v_matpower_new(2:end),'k-*'); hold off;

%% soc plot 
figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;

plot(0:number_of_market_period,[soc_init(1) soc_res(1,:)],'b-o','LineStyle', '--','DisplayName','central - BESS 1'); hold on;
plot(0:number_of_market_period,[soc_init(2) soc_res(2,:)],'c-o','LineStyle', '--','DisplayName','central - BESS 2'); hold on;
plot(0:number_of_market_period,[soc_init(3) soc_res(3,:)],'k-o','LineStyle', '--','DisplayName','central - BESS 3'); hold off;

legend('Location','NorthEast','Orientation','vertical');
legend boxoff 
legend show;

set(0,'DefaultAxesFontsSize',fntsize);
set(gcf,'Units','inche','Position',figSize);
xlabel('market period');
ylabel('SOC');

% filename='soc';
% export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');


%% consensus ADMM for LinDistFLow
%% clustering formation

% making clusters
mpc.cluster.number = 3;
%root cluster
mpc.cluster.root = 1; 
% cluster def
mpc.clusterBus{1} = [1:7 33:87 111]; % correct 
mpc.clusterBus{2} = [7 88:110]; 
mpc.clusterBus{3} = [7:32 112:141];


% root bus
mpc.rootBus{1} = 1;
mpc.rootBus{2} = 7;
mpc.rootBus{3} = 7;
% mpc.rootBus{4} = 15;
% mpc.rootBus{5} = 24;

mpc.CoupBus{1} = [7];
mpc.CoupBus{2} = [7];
mpc.CoupBus{3} = [7];


% PCC connection 
mpc.PCC{1} = 1;
mpc.PCC{2} = 0;
mpc.PCC{3} = 0;


% coulped bus index as slack bus 
mpc.CoupIdx{1} = 0;
mpc.CoupIdx{2} = 1;
mpc.CoupIdx{3} = 1;



%% network info extraction 
nbr = size(mpc.branch,1);
nbu = size(mpc.bus,1);
Ndg = size(mpc.gen,1);

% branch reformulation for cluster
% init flag
for c = 1:mpc.cluster.number
    mpc.cluster.branch_flag{c}=0;
    mpc.cluster.gen_flag{c}=0;
    mpc.cluster.strg_flag{c}=0;
end
% assign cluster branch & gennerations 
for c = 1:mpc.cluster.number
    %branch
    for br = 1:nbr
        if ismember(mpc.branch(br,F_BUS),mpc.clusterBus{c}) && ...
                ismember(mpc.branch(br,T_BUS),mpc.clusterBus{c})
            if mpc.cluster.branch_flag{c} == 0;            
                mpc.clusterBranch{c} = mpc.branch(br,:);
                mpc.cluster.branch_flag{c} = 1;
            else
                mpc.clusterBranch{c} = [mpc.clusterBranch{c}; mpc.branch(br,:)];
            end
        end 
    end
    
    % gennerations & gencost ( Price for P )
    for dg = 1:Ndg
        if ismember(mpc.gen(dg,1),mpc.clusterBus{c})
            if mpc.cluster.gen_flag{c} == 0;            
                mpc.clusterGen{c} = mpc.gen(dg,:);
                mpc.clusterGenCost{c} = mpc.gencost(dg,:);
                mpc.cluster.gen_flag{c} = 1;
            else
                mpc.clusterGen{c} = [mpc.clusterGen{c}; mpc.gen(dg,:)];
                mpc.clusterGenCost{c} = [mpc.clusterGenCost{c}; mpc.gencost(dg,:)]; 
            end
        end
    end
    
    % storage 
    for i_storage = 1:nstroage
        if ismember(mpc.storage(i_storage,1),mpc.clusterBus{c})
            if mpc.cluster.strg_flag{c} == 0;            
                mpc.clusterStorage{c} = mpc.storage(i_storage,:);

                mpc.cluster.strg_flag{c} = 1;
            else
                mpc.clusterStorage{c} = [mpc.clusterStorage{c}; mpc.storage(i_storage,:)];

            end
        end
    end
    
    
    % price for Q
    for dg = 1:Ndg
        if ismember(mpc.gen(dg,1),mpc.clusterBus{c})
                mpc.clusterGenCost{c} = [mpc.clusterGenCost{c}; mpc.gencost(Ndg+dg,:)]; % Price for P
        end
    end
end

% resign the branch index to internal idx(order same as in mpc.clusterBus)
for c = 1:mpc.cluster.number
    nbr_c = size(mpc.clusterBranch{c},1);
    for br = 1:nbr_c
        Idx_F_BUS = find(mpc.clusterBus{c}==mpc.clusterBranch{c}(br,F_BUS));
        mpc.clusterBranch{c}(br,F_BUS) = Idx_F_BUS;
        Idx_T_BUS = find(mpc.clusterBus{c}==mpc.clusterBranch{c}(br,T_BUS));
        mpc.clusterBranch{c}(br,T_BUS) = Idx_T_BUS;
    end
end

% resign the gen/storage index to internal idx (order same as in mpc.clusterBus)
for c = 1:mpc.cluster.number
    ndg_c = size(mpc.clusterGen{c},1);
    n_batt_c = size(mpc.clusterStorage{c},1);
    %gen
    for dg = 1:ndg_c
        Idx_GEN_BUS = find(mpc.clusterBus{c}==mpc.clusterGen{c}(dg,1));
        mpc.clusterGen{c}(dg,1) = Idx_GEN_BUS;
    end
    %storage
    for i_strg = 1:n_batt_c
        Idx_STRG_BUS = find(mpc.clusterBus{c}==mpc.clusterStorage{c}(i_strg,1));
        mpc.clusterStorage{c}(i_strg,1) = Idx_STRG_BUS;
    end
    
end

for c = 1:mpc.cluster.number
    nbr_c = size(mpc.clusterBranch{c},1);
    nbu_c = size(mpc.clusterBus{c},2);
    % save to mpc 
    mpc.cluster.nbu{c} = nbu_c;
    mpc.cluster.nbr{c} = nbr_c;       
end


% additional info for clusters 1
for c = 1:mpc.cluster.number
   mpc.cluster.InternalClusterBusIdx{c} = 1:mpc.cluster.nbu{c};
   mpc.cluster.genBus{c} = mpc.clusterGen{c}(:,1)';
end

%% mpc partitioning (for solving PF individually)
mpc_cluster_array = mpcPartition_v8_V2(mpc);

% new Ybus
for c = 1:mpc.cluster.number
    mpc.cluster.Ybus{c} = makeYbus(mpc_cluster_array(c));
    
    ReY_cluster = real(mpc.cluster.Ybus{c}(2:end,2:end));
    ImY_cluster = imag(mpc.cluster.Ybus{c}(2:end,2:end));
    
    mpc.cluster.LinDist_c{c} = inv([ReY_cluster -ImY_cluster; -ImY_cluster -ReY_cluster]);
    
    %% construct the linearized line current model
    
    mpc.cluster.Y_line_abs{c} = zeros(mpc.cluster.nbr{c},mpc.cluster.nbr{c});    
    mpc.cluster.Cf{c} = zeros(mpc.cluster.nbr{c},mpc.cluster.nbu{c});
    mpc.cluster.Ct{c} = zeros(mpc.cluster.nbr{c},mpc.cluster.nbu{c}); 
    
    for br = 1:mpc.cluster.nbr{c}
        br_F_BUS = mpc_cluster_array(c).branch(br,F_BUS);
        br_T_BUS = mpc_cluster_array(c).branch(br,T_BUS);        
        br_BR_R = mpc_cluster_array(c).branch(br,BR_R);
        br_BR_X = mpc_cluster_array(c).branch(br,BR_X);
        br_Y = 1 / (br_BR_R + 1j * br_BR_X);    
        
        mpc.cluster.Y_line_abs{c}(br,br) = abs(br_Y);       
        mpc.cluster.Cf{c}(br,br_F_BUS) = 1;
        mpc.cluster.Ct{c}(br,br_T_BUS) = 1;
    end
        
end

tic; 
%% ADMM init 
ADMM.parameter.num_itr = 200;
ADMM.parameter.mu = 10;
ADMM.parameter.tau = 0.2;

% ADMM varibales declaration 
%ADMM.var.global_bus_set = [7 15 24];
Ngc = 1;

%% % third dimension is for the time 
for t = 1:time_step  
    ADMM.var.p_global{t} = zeros(mpc.cluster.number,Ngc);
    ADMM.var.q_global{t} = zeros(mpc.cluster.number,Ngc);
    ADMM.var.v_global{t} = zeros(mpc.cluster.number,Ngc);

    % previous value for dual residual calucalation 
    ADMM.var.p_global_pre{t} = zeros(1,Ngc);
    ADMM.var.q_global_pre{t} = zeros(1,Ngc);
    ADMM.var.v_global_pre{t} = zeros(1,Ngc);
    
    % temp variable for averaging step 
    ADMM.var.v_global_temp{t} = zeros(mpc.cluster.number,Ngc);
    ADMM.var.p_global_temp{t} = zeros(mpc.cluster.number,Ngc);
    ADMM.var.q_global_temp{t} = zeros(mpc.cluster.number,Ngc); 
    
    lambda_rec{t} = zeros(3,ADMM.parameter.num_itr);
    kappa_rec{t} = zeros(3,ADMM.parameter.num_itr);
%     P23{t} = [];
%     Q23{t} = [];
%     P13{t} = [];
%     Q13{t} = [];
%     P12{t} = [];
%     Q12{t} = [];
end
       
% dimension of Lambda - cluster 1: 2xtime_step, cluster 2: 2xtime_step, cluster 3: 4xtime_step,
% cluster 4: 4xtime_step, cluster 5:2xtime_step
ADMM.var.Lambda = cell(mpc.cluster.number,time_step);

ADMM.var.Lambda_pre = cell(mpc.cluster.number,time_step);

ADMM.var.Voltage = cell(mpc.cluster.number,time_step); % reference voltage for each cluster
ADMM.var.VoltageProfile = cell(mpc.cluster.number,1); % nodal voltage vector of each cluster
ADMM.var.p_gen = cell(mpc.cluster.number,time_step);
ADMM.var.q_gen= cell(mpc.cluster.number,time_step);
%ADMM.var.Lambda_va = cell(mpc.cluster.number,time_step); % voltage and angles are exact the same for PQ injections 
ADMM.norm2aux = zeros(1,ADMM.parameter.num_itr);
ADMM.var.rho = 8e-2*ones(mpc.cluster.number,1);%8e-2 for t = 24 0.05
rho_init = ADMM.var.rho;

%primal residual declare
ADMM.var.r = zeros(mpc.cluster.number,1);
% dual residual declare
ADMM.var.s= zeros(mpc.cluster.number,1);

ADMM.var.SOC = cell(mpc.cluster.number,1);

for c = 1:mpc.cluster.number
        NofCoupBus = size(mpc.CoupBus{c},2);  
        % lambda init same order as in CoupBus
        for t = 1:time_step            
            if c == 1
                ADMM.var.Lambda{c,t} = 0e1*ones(NofCoupBus,2);
                ADMM.var.Voltage{c,t} = mpc_cluster_array(c).gen(1,VG); 
            else 
                ADMM.var.Lambda{c,t} = 0e1*ones(NofCoupBus,2);
                ADMM.var.Voltage{c,t} = mpc_cluster_array(c).gen(1,VG);
            end         
        end 
end

% % local objective values 
% obj_v = zeros(1,ADMM.parameter.num_itr);
obj_v_temp = zeros(3,1);

% fast ADMM 
alpha = ones(1,time_step);
ADMM.var.Lambda_pre = ADMM.var.Lambda;
ADMM.var.p_global_histroy = cell(ADMM.parameter.num_itr,time_step);
ADMM.var.q_global_histroy = cell(ADMM.parameter.num_itr,time_step);

%%
yalmip('clear');

%cost_debug = zeros(5,4);

for i = 1:ADMM.parameter.num_itr   	
    %% assign new operational point 
    for c = 1:mpc.cluster.number
            sprintf('ADMM : %d of %d interations, cluster: %d', [i,ADMM.parameter.num_itr, c])    
            %% new loop 
            soc_init_cluster = soc_init(c);

            Bess_capacity = mpc_cluster_array(c).storage(:,2);
            Bess_charge_efficiency = mpc_cluster_array(c).storage(:,3);
            Bess_discharge_efficiency = mpc_cluster_array(c).storage(:,4);
            Bess_cost_coefficient = mpc_cluster_array(c).storage(:,6);
          
            % injection constr
            QG_max = zeros(mpc.cluster.nbu{c},1);
            PG_max = zeros(mpc.cluster.nbu{c},1);
            QG_min = zeros(mpc.cluster.nbu{c},1);
            PG_min = zeros(mpc.cluster.nbu{c},1);
            
            QG_max(mpc_cluster_array(c).gen(:,1)) = mpc_cluster_array(c).gen(:,QMAX)/mpc.baseMVA;
            PG_max(mpc_cluster_array(c).gen(:,1)) = mpc_cluster_array(c).gen(:,PMAX)/mpc.baseMVA;
            QG_min(mpc_cluster_array(c).gen(:,1)) = mpc_cluster_array(c).gen(:,QMIN)/mpc.baseMVA;
            PG_min(mpc_cluster_array(c).gen(:,1)) = mpc_cluster_array(c).gen(:,PMIN)/mpc.baseMVA;
            
            % ramping rate limit of energy storage 
            PG_max(mpc_cluster_array(c).storage(:,1)) = mpc_cluster_array(c).storage(:,5)/mpc.baseMVA;
            PG_min(mpc_cluster_array(c).storage(:,1)) = -mpc_cluster_array(c).storage(:,5)/mpc.baseMVA;
            
            % SOC limit 
%             SOC_min = 0.3;
%             SOC_max = 0.99;
            
            % Generation Set points
            PG_var = sdpvar(mpc.cluster.nbu{c},time_step);
            QG_var = sdpvar(mpc.cluster.nbu{c},time_step);
            
            % only one battery allowed
            P_charge_var = sdpvar(1,time_step);
            P_discharge_var = sdpvar(1,time_step);
            SOC_var = sdpvar(1,time_step);

            PQ_inj_var = [PG_var(2:end,:) - repmat(mpc_cluster_array(c).bus(2:end,PD)/mpc.baseMVA,1,time_step)*diag(load_profile); QG_var(2:end,:) - repmat(mpc_cluster_array(c).bus(2:end,QD)/mpc.baseMVA,1,time_step)*diag(load_profile)];
            deltaVA = - mpc.cluster.LinDist_c{c}*PQ_inj_var; % assume the load does not change 
                                              
            PQ0_constr = [];
            V_constr = [];
            PQ_constr = [];
            Batt_constr = [];
            SOC_constr = [];

            P0 = sum(mpc_cluster_array(c).bus(:,PD))/mpc.baseMVA;
            Q0 = sum(mpc_cluster_array(c).bus(:,QD))/mpc.baseMVA;
            
%             if c~= 1
%                 Delta_v0 = sdpvar(1,time_step); % p.u             % v0 as variable 
%             end
            
            vbase = zeros(mpc.cluster.nbu{c}-1,time_step);
            for t = 1:time_step
                % base voltage 
%                 if c~= 1
%                     vbase(:,t) = (mpc_cluster_array(c).gen(1,VG)+Delta_v0(t))*ones(mpc.cluster.nbu{c}-1,1);
%                 else
                    vbase(:,t) = ADMM.var.Voltage{c,t}*ones(mpc.cluster.nbu{c}-1,1);
%                 end
                
            end
            % 
            VM_var = deltaVA(1:(mpc.cluster.nbu{c}-1),:) + vbase;
            % constraint time loop
            for t = 1:time_step
                % p_0, q_0 constr
                PQ0_constr = [PQ0_constr, PG_var(1,t) == P0*load_profile(t)-sum(PG_var(2:end,t)), QG_var(1,t) == Q0*load_profile(t)-sum(QG_var(2:end,t))];
                % voltage constraints        
%                 if c~= 1
%                     V_constr = [V_constr, 0.95 <= mpc_cluster_array(c).gen(1,VG)+Delta_v0(t) <= 1.05 ,0.95*ones(mpc.cluster.nbu{c}-1,1)/1.001 <= VM_var(:,t) <= 1.05*ones(mpc.cluster.nbu{c}-1,1)*1.00]; % 5% maximal error is tolerated
%                 else
                V_constr = [V_constr, 0.90*ones(mpc.cluster.nbu{c}-1,1)/1.001 <= VM_var(:,t) <= 1.05*ones(mpc.cluster.nbu{c}-1,1)*1]; % 5% maximal error is tolerated % it doesn't mater
%                 end
                %%% stopped here
                      
                PQ_constr = [PQ_constr, PG_min(2:end)<=PG_var(2:end,t)<=PG_max(2:end),QG_min(2:end)<=QG_var(2:end,t)<=QG_max(2:end), PG_min(1)<= PG_var(1,t) <= PG_max(1), QG_min(1)<= QG_var(1,t) <= QG_max(1)];
                Batt_constr = [Batt_constr, PG_var(mpc_cluster_array(c).storage(:,1),t) == P_discharge_var(t) - P_charge_var(t), P_discharge_var(t)>=0, P_charge_var(t) >= 0];

                if t == 1
                    SOC_constr = [SOC_var(t) == soc_init_cluster - inv(diag(Bess_capacity))*(inv(diag(Bess_discharge_efficiency))*P_discharge_var(t) - inv(diag(Bess_charge_efficiency))*P_charge_var(t)), SOC_var(t)>=SOC_min, SOC_var(t)<=SOC_max];
                else
                    SOC_constr = [SOC_constr, SOC_var(t) == SOC_var(t-1) - inv(diag(Bess_capacity))*(inv(diag(Bess_discharge_efficiency))*P_discharge_var(t) - inv(diag(Bess_charge_efficiency))*P_charge_var(t)), SOC_var(t)>=SOC_min, SOC_var(t)<=SOC_max];
                end
            end
            
            %% robust optimization varibales 
            mu_sum = sdpvar(1,1);
            mu_prox = sdpvar(1,time_step);
            % robust constraint 
            robust_constr = [];
            % robustification constraint 1
            robust_constr = [mu_sum>=0];

            % robustification constraint 2
            for t = 1:time_step
                robust_constr = [robust_constr, mu_sum + mu_prox(t) >= max_deviation_from_psp_price(t)*(PG_var(1,t)), mu_prox(t)>= 0]; 
            end

            
            if c == 1 
                constr = [PQ0_constr,V_constr,PQ_constr,Batt_constr,SOC_constr,robust_constr];
            else 
                constr = [PQ0_constr,V_constr,PQ_constr,Batt_constr,SOC_constr];
            end
           
%            congestion constraints 
%             % line current constr
%            if c~= 1
%                  v_all_bus_var = [v0_var; deltaVA(1:(mpc.cluster.nbu{c}-1)) + vbase]; % 5% maximal error is tolerated
%            else
%                  v_all_bus_var = [mpc_cluster_array(c).gen(1,VG); deltaVA(1:(mpc.cluster.nbu{c}-1)) + vbase]; % 5% maximal error is tolerated
%            end
% 
%            line_current_magnitude = mpc.cluster.Y_line_abs{c}*(mpc.cluster.Cf{c}*v_all_bus_var - mpc.cluster.Ct{c}*v_all_bus_var);
    
%            if  c == 3
%                congestion_constr = [-0.05*sqrt(2)<line_current_magnitude(1)<=0.05*sqrt(2)];
%                Constr = [PQ0_constr,V_constr,PQ_constr,congestion_constr];
%            else
%                Constr = [PQ0_constr,V_constr,PQ_constr];
%            end
% 
%             % voltage regulation in the objective
%             v_norm = ones(mpc.cluster.nbu{c}-1,1);
%             weighting_factor = 5*ones(mpc.cluster.nbu{c}-1,1);
%             v_regulation = (deltaVA(1:(mpc.cluster.nbu{c}-1)) + vbase - v_norm)'*diag(weighting_factor)*(deltaVA(1:(mpc.cluster.nbu{c}-1)) + vbase - v_norm);
            
            % objective    
%             dg_cost = 0;
%             battery_cost = 0;
%             psp_cost = 0;
%             
%             for t = 1:time_step
%                 dg_cost = dg_cost + (mpc.gencost(2:ngen,5)'*PG_var(mpc.gen(2:end,1),t) + mpc.gencost(ngen+2:end,5)'*QG_var(mpc.gen(2:end,1),t)); 
%                 battery_cost = battery_cost + Bess_cost_coefficient'*(P_discharge_var(:,t)-P_charge_var(:,t)); 
%                 psp_cost = psp_cost + (PG_var(1,t)+P0)*psp_price(t) + (QG_var(1,t)+Q0)*psp_price(t)/3;
%             end
% 
%             cost = battery_cost + psp_cost/7 + dg_cost; % 7 to scale down the price to normal price

        robust_cost = 0;
        coupled_cost = 0;
        dg_cost = 0;
        battery_cost = 0;
        psp_cost = 0;      
        total_cost = 0;  
        
        
        if c == 1 % cluster 1 yes
                % robustification objective 
                for t = 1:time_step
                    robust_cost = robust_cost + (mu_prox(t) + (PG_var(1,t))*psp_price(t));
                end

                robust_cost = robust_cost + level_of_conserv*mu_sum;
            
                for t = 1:time_step
                    coupled_cost = coupled_cost + ADMM.var.Lambda{c,t}(1)*(PG_var(7,t) - (- ADMM.var.p_global{t}(c,1))) + ...  % 
                    0.5*ADMM.var.rho(c)*(PG_var(7,t) - (-ADMM.var.p_global{t}(c,1)))^2 + ...
                    ADMM.var.Lambda{c,t}(2)*(QG_var(7,t) - (-ADMM.var.q_global{t}(c,1))) + ...
                    0.5*ADMM.var.rho(c)*(QG_var(7,t) - (-ADMM.var.q_global{t}(c,1)))^2;
                
                    dg_cost = dg_cost + (mpc_cluster_array(c).gencost(2:3,5)'*PG_var(mpc_cluster_array(c).gen(2:3,1),t)+ PG_var(mpc_cluster_array(c).gen(2:3,1),t)'*diag(mpc_cluster_array(c).gencost(2:3,6))'*PG_var(mpc_cluster_array(c).gen(2:3,1),t) + mpc_cluster_array(c).gencost(6:7,5)'*QG_var(mpc_cluster_array(c).gen(2:3,1),t)); 
                    battery_cost = battery_cost + Bess_cost_coefficient'*(P_discharge_var(t)-P_charge_var(t))^2; 
                    psp_cost = psp_cost + (PG_var(1,t))*psp_price(t) + (QG_var(1,t))*psp_price(t)/3;                
                end
                
                total_cost = coupled_cost + (dg_cost*price_factor + battery_cost*price_factor + psp_cost) + robust_cost;
                
        elseif c == 2   % cluster 2 & 5 yes

                %% stopped here
                for t = 1:time_step
                    coupled_cost = coupled_cost + ...
                    ADMM.var.Lambda{c,t}(1)*(PG_var(1,t) - ADMM.var.p_global{t}(c,1)) + ...% 
                    0.5*ADMM.var.rho(c)*(PG_var(1,t) - ADMM.var.p_global{t}(c,1))^2 + ...% 
                    ADMM.var.Lambda{c,t}(2)*(QG_var(1,t) - ADMM.var.q_global{t}(c,1)) + ...% 
                    0.5*ADMM.var.rho(c)*(QG_var(1,t) - ADMM.var.q_global{t}(c,1))^2;
               
                    dg_cost = dg_cost + (mpc_cluster_array(c).gencost(2,5)'*PG_var(mpc_cluster_array(c).gen(2:end,1),t) + mpc_cluster_array(c).gencost(4,5)'*QG_var(mpc_cluster_array(c).gen(2:end,1),t));
                    
                    battery_cost = battery_cost + Bess_cost_coefficient'*(P_discharge_var(t)-P_charge_var(t))^2; 
                    
                end
                
                total_cost = coupled_cost + (dg_cost*price_factor + battery_cost*price_factor);
                
        elseif c == 3 %  cluster 3 yes 
            
            for t = 1:time_step  
                % lambda order [c1:p, c1:q; c2:p, c2:q]
                coupled_cost = coupled_cost + ...
                ADMM.var.Lambda{c,t}(1,1)*(PG_var(1,t) - ADMM.var.p_global{t}(c,1)) + ...% 
                0.5*ADMM.var.rho(c)*(PG_var(1,t) - ADMM.var.p_global{t}(c,1))^2 + ...%
                ADMM.var.Lambda{c,t}(1,2)*(QG_var(1,t) - ADMM.var.q_global{t}(c,1)) + ...%
                0.5*ADMM.var.rho(c)*(QG_var(1,t) - ADMM.var.q_global{t}(c,1))^2;

            
                dg_cost = dg_cost + (mpc_cluster_array(c).gencost(2:4,5)'*PG_var(mpc_cluster_array(c).gen(2:end,1),t) + PG_var(mpc_cluster_array(c).gen(2:end,1),t)'*diag(mpc_cluster_array(c).gencost(2:4,6))'*PG_var(mpc_cluster_array(c).gen(2:end,1),t) + mpc_cluster_array(c).gencost(6:end,5)'*QG_var(mpc_cluster_array(c).gen(2:end,1),t)); 
                battery_cost = battery_cost + Bess_cost_coefficient'*(P_discharge_var(t)-P_charge_var(t))^2;            
            end
            
            total_cost = coupled_cost + (dg_cost*price_factor + battery_cost*price_factor);

        end

        %Obj = Obj + v_regulation;
        %ops = sdpsettings('solver','gurobi');
        optimize(constr,total_cost,ops);

        %% results extraction        
        Pvar_res = value(PG_var);
        PGres = Pvar_res(mpc_cluster_array(c).gen(:,1),:); % EXTRACT THE NECESSARY PART 
        Qvar_res = value(QG_var);
        QGres = Qvar_res(mpc_cluster_array(c).gen(:,1),:);
        Vm_res = value(VM_var);        
        v_lin_new = Vm_res;
        
        obj_v_temp(c) = value(total_cost);
        
%         %% new linDist Power Flow
%         u_lin_new = mpc.cluster.LinDist_c{c}*[-Pvar_res(2:end) + mpc_cluster_array(c).bus(2:end,PD)/mpc.baseMVA;  -Qvar_res(2:end) + mpc_cluster_array(c).bus(2:end,QD)/mpc.baseMVA];
%         
%         if c~= 1
%             v_lin_new = -u_lin_new(1:(mpc.cluster.nbu{c}-1)) + value(vbase);
%             v0_res = value(v0_var);
%         else
%             v_lin_new = -u_lin_new(1:(mpc.cluster.nbu{c}-1)) + vbase;  
%         end
 
%        ADMM.var.Voltage{c} = v_lin_new;
        
        ADMM.var.VoltageProfile{c} = Vm_res;
        ADMM.var.p_gen{c} = PGres;
        ADMM.var.q_gen{c} = QGres;
        ADMM.var.SOC{c} = value(SOC_var);
        
        
        cost_debug(c,1) = price_factor*value(battery_cost);
        cost_debug(c,2) = value(psp_cost);
        cost_debug(c,3) = price_factor*value(dg_cost);
        cost_debug(c,4) = value(coupled_cost);
%         
        %%%%% price decomposition 
%         if i == ADMM.parameter.num_itr
%             if c == 1
%                 Vbase = [1; v_lin_new];
%                 %% debug 
%                 %test = runpf(mpc_cluster_array(1));
%                 %Abase = test.bus(1:end,VA)/180*pi;
%                 %%
%                 Abase = [0; -u_lin_new(mpc.cluster.nbu{c}:end)]; % in radius 
%                 Ybus =  mpc.cluster.Ybus{c};
%                 n = size(mpc_cluster_array(c).bus,1);
%                 nnb = n-1;
%                                 
%                 % Ploss0 and Qloss0 are in per unit (can be compared to sum(get_losses(PFres)))
%                 [Agrid2, Aloss, Ploss0, Qloss0,] = getLinModel_clean(Vbase, Abase, Ybus, n, mpc_cluster_array(c));
%                 [LF_lin_dPloss_dPd, LF_lin_dPloss_dQd,LF_lin_dQloss_dPd, LF_lin_dQloss_dQd, Jlin] = getLF(mpc_cluster_array(c), n, Agrid2, Aloss) ;
%                 [LF_lin_dPloss_dPd, LF_lin_dPloss_dQd,LF_lin_dQloss_dPd, LF_lin_dQloss_dQd, Jlin, Jred] = getLF_v2(mpc_cluster_array(c), n, Agrid2, Aloss); 
%              
%                 mu_v = dual(V_constr);   
%                 mu_v_max = mu_v(1:nnb);
%                 mu_v_min = mu_v(nnb+1:end);
% 
%                 % Get DLMPs
%                 Kvp = mpc.cluster.LinDist_c{c}(1:nnb,1:nnb);
%                 Kvq = mpc.cluster.LinDist_c{c}(1:nnb,n:end);
%                                 
%                 Mu_V_diff =  mu_v_max - mu_v_min;
%                 DLMP_V = Kvp'*Mu_V_diff;
%                 DLMP_V_q = Kvq'*Mu_V_diff;
%                                            
%                 % P Q DLMPs calculation 
%                 P_DLMPs = 17 - (17)*LF_lin_dPloss_dPd(2:end) - 3*LF_lin_dQloss_dPd(2:end) - DLMP_V/mpc.baseMVA;
%                 Q_DLMPs = 3 - (17)*LF_lin_dPloss_dQd(2:end) - 3*LF_lin_dQloss_dQd(2:end) - DLMP_V_q/mpc.baseMVA;
%             end
%             
%                if  c == 3
%                    congestion_dual = dual(congestion_constr);
%                    
%                    congestion_dual_diff_br1 = congestion_dual(2)-congestion_dual(1);
%                    congestion_dual_diff_full = zeros(mpc.cluster.nbr{3},1);
%                    congestion_dual_diff_full(1) = congestion_dual_diff_br1;
%                                     
%                    KKK = mpc.cluster.LinDist_c{3};
%                    Kvp = KKK(1:(mpc.cluster.nbu{3}-1),1:(mpc.cluster.nbu{3}-1));
%                    Kvq = KKK(1:(mpc.cluster.nbu{3}-1),mpc.cluster.nbu{3}:end);
%                    
%                    
%                    temp = (mpc.cluster.Cf{3} - mpc.cluster.Ct{3})'*mpc.cluster.Y_line_abs{3}'*congestion_dual_diff_full/10;
%                    
%                    congestion_part_trading_price_p = -Kvp*temp(2:end);
%                    congestion_part_trading_price_q = -Kvq*temp(2:end);
%                                                       
%                end
%             
%         end
        %% local ADMM update - save the results for trust region     
        for t = 1:time_step
            if c == 1 % upstream 
                ADMM.var.v_global_temp{t}(c,1) = v_lin_new(6,t); %for averaging purpose 

                ADMM.var.p_global_temp{t}(c,1) = -PGres(4,t); %ADMM update 
                ADMM.var.q_global_temp{t}(c,1) = -QGres(4,t); %ADMM update 
                
            elseif c== 2 

                ADMM.var.p_global_temp{t}(c,mpc.CoupIdx{c}) = PGres(1,t) ; % from new PF point % p.u 
                ADMM.var.q_global_temp{t}(c,mpc.CoupIdx{c}) = QGres(1,t) ;

%                ADMM.var.v_global_temp(c,temp) = v0_res; 
            elseif c == 3
                %coupled bus 1
                ADMM.var.p_global_temp{t}(c,1) = PGres(1,t); 
                ADMM.var.q_global_temp{t}(c,1) = QGres(1,t);    
            end
        end
        
            %% clear memory    
            yalmip('clear');
    end   
                      

         
    %% global ADMM update 
    % Averaging step manuelly....
    % Area 1
    for t = 1:time_step
        P23 = sum(ADMM.var.p_global_temp{t}(2:3,1));
        temp1 = 0.5*(P23 + ADMM.var.p_global_temp{t}(1,1));    
        Q23 = sum(ADMM.var.q_global_temp{t}(2:3,1));
        temp4 = 0.5*(Q23 + ADMM.var.q_global_temp{t}(1,1));

        % Area 2 
        P13 = ADMM.var.p_global_temp{t}(1,1) - ADMM.var.p_global_temp{t}(3,1);
        temp2 = 0.5*(P13 + ADMM.var.p_global_temp{t}(2,1));    
        Q13 = ADMM.var.q_global_temp{t}(1,1) - ADMM.var.q_global_temp{t}(3,1);
        temp5 = 0.5*(Q13 + ADMM.var.q_global_temp{t}(2,1));

        % Area 3 
        P12 = ADMM.var.p_global_temp{t}(1,1) - ADMM.var.p_global_temp{t}(2,1);
        temp3 = 0.5*(P12 + ADMM.var.p_global_temp{t}(3,1));   
        Q12 = ADMM.var.q_global_temp{t}(1,1) - ADMM.var.q_global_temp{t}(2,1);
        temp6 = 0.5*(Q12 + ADMM.var.q_global_temp{t}(3,1));


        % assign p & q coupled bus 1
        ADMM.var.p_global{t}(1,1) = temp1; % in p.u 
        ADMM.var.p_global{t}(2,1) = temp2; % in p.u 
        ADMM.var.p_global{t}(3,1) = temp3; % in p.u       
        ADMM.var.q_global{t}(1,1) = temp4; % in p.u 
        ADMM.var.q_global{t}(2,1) = temp5; % in p.u 
        ADMM.var.q_global{t}(3,1) = temp6; % in p.u   

        
        % update reference voltage 
        ADMM.var.Voltage{2,t} = ADMM.var.v_global_temp{t}(1,1);
        ADMM.var.Voltage{3,t} = ADMM.var.v_global_temp{t}(1,1);
        
        
        % for fast ADMM
        ADMM.var.p_global_histroy{i,t} = ADMM.var.p_global{t};
        ADMM.var.q_global_histroy{i,t} = ADMM.var.q_global{t};
         
    end
                      
%     %% updating primal/dual residuals % not including voltage yet
    if i > 1
        diff_p_tempPrim = 0;
        diff_q_tempPrim = 0;
        diff_p_tempDual = 0;
        diff_q_tempDual = 0;
        
        temp = 0;
        
        for t = 1:time_step
            diff_p_tempPrim = diff_p_tempPrim + sum((ADMM.var.p_global{t} -  ADMM.var.p_global_temp{t}).^2,2);
            diff_q_tempPrim = diff_q_tempPrim + sum((ADMM.var.q_global{t} -  ADMM.var.q_global_temp{t}).^2,2);
    %        diff_v_tempPrim = sum((ADMM.var.v_global -  ADMM.var.v_global_temp).^2,2);
            %         
    %         % record the norm of disagreement
            %ADMM.norm2aux(:,i) = [diff_p_tempPrim;diff_q_tempPrim];
    %         
            diff_p_tempDual = diff_p_tempDual + sum((ADMM.var.p_global{t} -  ADMM.var.p_global_pre{t}).^2,2);
            diff_q_tempDual = diff_q_tempDual + sum((ADMM.var.q_global{t} -  ADMM.var.q_global_pre{t}).^2,2);
    %        diff_v_tempDual = sum((ADMM.var.v_global -  ADMM.var.v_global_pre).^2,2);
    
        end
        
             ADMM.var.r = diff_p_tempPrim + diff_q_tempPrim; % prime residual square
             ADMM.norm2aux(i) = sum(ADMM.var.r);
             ADMM.var.s = ADMM.var.rho.^2.*(diff_p_tempDual + diff_q_tempDual);  % dual residual square
             
    end
    
%     if i > 1
%         ADMM.var.Lambda_temp_previous = ADMM.var.Lambda_temp;
%     end 
    
%     %% fast ADMM 
%     % lambda update 
%     for t = 1:time_step
%         ADMM.var.Lambda_temp{1,t}(1) = ADMM.var.Lambda{1,t}(1) - ADMM.var.rho(1)*(ADMM.var.p_global_temp{t}(1,1) - ADMM.var.p_global{t}(1,1));
%         ADMM.var.Lambda_temp{1,t}(2) = ADMM.var.Lambda{1,t}(2) - ADMM.var.rho(1)*(ADMM.var.q_global_temp{t}(1,1) - ADMM.var.q_global{t}(1,1));
% 
%         % downstream 
%         c = 2; 
%         ADMM.var.Lambda_temp{c,t}(1) = ADMM.var.Lambda{c,t}(1) + ADMM.var.rho(c)*(ADMM.var.p_global_temp{t}(2,1) - ADMM.var.p_global{t}(c,1));
%         ADMM.var.Lambda_temp{c,t}(2) = ADMM.var.Lambda{c,t}(2) + ADMM.var.rho(c)*(ADMM.var.q_global_temp{t}(2,1) - ADMM.var.q_global{t}(c,1));
%         
%         c = 3;
%         ADMM.var.Lambda_temp{c,t}(1,1) = ADMM.var.Lambda{c,t}(1,1) + ADMM.var.rho(c)*(ADMM.var.p_global_temp{t}(c,1) - ADMM.var.p_global{t}(c,1));
%         ADMM.var.Lambda_temp{c,t}(1,2) = ADMM.var.Lambda{c,t}(1,2) + ADMM.var.rho(c)*(ADMM.var.q_global_temp{t}(c,1) - ADMM.var.q_global{t}(c,1));
% 
%     end
    
    
  %  alpha(i+1) = 0.5*(1 + sqrt(1+4*(alpha(i))^2));
    
    % update previous global var 
    ADMM.var.p_global_pre = ADMM.var.p_global;
    ADMM.var.q_global_pre = ADMM.var.q_global;
%    ADMM.var.v_global_pre = ADMM.var.v_global;
%     if i > 1
%         for t = 1:time_step
%             for c = 1:3
%                 if c == 1
%                     ADMM.var.Lambda{c,t}  =  ADMM.var.Lambda_temp{c,t} -  alpha(i-1)/alpha(i+1)*(ADMM.var.Lambda_temp{c,t} -  ADMM.var.Lambda_temp_previous{c,t});     
%                 else 
%                     ADMM.var.Lambda{c,t}  =  ADMM.var.Lambda_temp{c,t} +   alpha(i-1)/alpha(i+1)*(ADMM.var.Lambda_temp{c,t} -  ADMM.var.Lambda_temp_previous{c,t}); 
%                 end
%             end 
%         end
%         
%         for t = 1:time_step
%                 ADMM.var.p_global{t}  =  ADMM.var.p_global{t} +  alpha(i-1)/alpha(i+1)*(ADMM.var.p_global_histroy{i,t} - ADMM.var.p_global_histroy{i-1,t});     %% two cases if i = 1
%                 ADMM.var.q_global{t}  =  ADMM.var.q_global{t} +  alpha(i-1)/alpha(i+1)*(ADMM.var.q_global_histroy{i,t} - ADMM.var.q_global_histroy{i-1,t}); 
%         end
%               
%     else 
%         for t = 1:time_step
%             for c = 1:3
%                 ADMM.var.Lambda{c,t}  =  ADMM.var.Lambda_temp{c,t};     %% p 
%             end
%         end       
%     end
    
   if i <= ADMM.parameter.num_itr && i > 1 % this is the correct code here, but it works much worse than the wrong one!
        for c = 1:mpc.cluster.number             
            if ADMM.var.r(c)>ADMM.parameter.mu*ADMM.var.s(c)
                ADMM.var.rho(c) = ADMM.var.rho(c)*(1+ADMM.parameter.tau);
            elseif ADMM.var.s(c)>ADMM.parameter.mu*ADMM.var.r(c)
                ADMM.var.rho(c) = ADMM.var.rho(c)/(1+ADMM.parameter.tau);
            else
                ADMM.var.rho(c) = ADMM.var.rho(c);
            end
        end
   else  
       %ADMM.var.rho = rho_init;
   end    
% %     


    %% multiplier updating step
    % upstream c=1
    for t = 1:time_step
        ADMM.var.Lambda{1,t}(1) = ADMM.var.Lambda{1,t}(1) - ADMM.var.rho(1)*(ADMM.var.p_global_temp{t}(1,1) - ADMM.var.p_global{t}(1,1));
        ADMM.var.Lambda{1,t}(2) = ADMM.var.Lambda{1,t}(2) - ADMM.var.rho(1)*(ADMM.var.q_global_temp{t}(1,1) - ADMM.var.q_global{t}(1,1));

    %      ADMM.var.Lambda_va{1}(1) = ADMM.var.Lambda_va{1}(1) + ADMM.var.rho(1)*(ADMM.var.v_global_temp(1,1) - ADMM.var.v_global(1,1));

        % downstream 
        c = 2; 
        ADMM.var.Lambda{c,t}(1) = ADMM.var.Lambda{c,t}(1) + ADMM.var.rho(c)*(ADMM.var.p_global_temp{t}(2,1) - ADMM.var.p_global{t}(c,1));
        ADMM.var.Lambda{c,t}(2) = ADMM.var.Lambda{c,t}(2) + ADMM.var.rho(c)*(ADMM.var.q_global_temp{t}(2,1) - ADMM.var.q_global{t}(c,1));

    %     ADMM.var.Lambda_va{c}(1) = ADMM.var.Lambda_va{c}(1) + ADMM.var.rho(c)*(ADMM.var.v_global_temp(c,1) - ADMM.var.v_global(c,1));

        c = 3;
        ADMM.var.Lambda{c,t}(1,1) = ADMM.var.Lambda{c,t}(1,1) + ADMM.var.rho(c)*(ADMM.var.p_global_temp{t}(c,1) - ADMM.var.p_global{t}(c,1));
        ADMM.var.Lambda{c,t}(1,2) = ADMM.var.Lambda{c,t}(1,2) + ADMM.var.rho(c)*(ADMM.var.q_global_temp{t}(c,1) - ADMM.var.q_global{t}(c,1));


    end
    
%     ADMM.var.Lambda_va{c}(1) = ADMM.var.Lambda_va{c}(1) + ADMM.var.rho(c)*(ADMM.var.v_global_temp(c,3) - ADMM.var.v_global(c,3));
    
    %recording value of generations 
    %PQG_rec(:,i) = [TRopt.storage.PFsolution{1}.gen(1:end-1,PG); TRopt.storage.PFsolution{2}.gen(2:end,PG);TRopt.storage.PFsolution{3}.gen(2:end-1,PG);TRopt.storage.PFsolution{4}.gen(2:end-1,PG);TRopt.storage.PFsolution{5}.gen(2:end,PG);...
        %TRopt.storage.PFsolution{1}.gen(1:end-1,QG); TRopt.storage.PFsolution{2}.gen(2:end,QG);TRopt.storage.PFsolution{3}.gen(2:end-1,QG);TRopt.storage.PFsolution{4}.gen(2:end-1,QG);TRopt.storage.PFsolution{5}.gen(2:end,QG);];
   for t = 1:time_step 

        lambda_rec{t}(1,i) = ADMM.var.Lambda{1,t}(1);
        lambda_rec{t}(2,i) = ADMM.var.Lambda{2,t}(1);
        lambda_rec{t}(3,i) = ADMM.var.Lambda{3,t}(1,1);
        
        kappa_rec{t}(1,i) = ADMM.var.Lambda{1,t}(2);
        kappa_rec{t}(2,i) = ADMM.var.Lambda{2,t}(2);
        kappa_rec{t}(3,i) = ADMM.var.Lambda{3,t}(1,2);

             
   end
% %     sigma_rec(1,i) = ADMM.var.Lambda_va{1}(1);
% %     sigma_rec(2,i) = ADMM.var.Lambda_va{2}(1);
% %     sigma_rec(3,i) = ADMM.var.Lambda_va{3}(1,1);
% %     sigma_rec(4,i) = ADMM.var.Lambda_va{3}(2,1);
% %     sigma_rec(5,i) = ADMM.var.Lambda_va{4}(1,1);
% %     sigma_rec(6,i) = ADMM.var.Lambda_va{4}(2,1);
% %     sigma_rec(7,i) = ADMM.var.Lambda_va{5}(1);
% %     
    obj_v(i) = sum(obj_v_temp);   
%     mpc_cluster_array(2).gen(1,VG) = 1;
%     mpc_cluster_array(3).gen(1,VG) = 1;
%     mpc_cluster_array(4).gen(1,VG) = 1;
%     mpc_cluster_array(5).gen(1,VG) = 1;
    %mu_admm(i) = mu_v;   
end
%% plot energy inport/export 
p_inout = zeros(7,time_step);
q_inout = zeros(7,time_step);
for t = 1:time_step 
    p_inout(1,t) =  -mpc.baseMVA*ADMM.var.p_global{t}(1,1);
    p_inout(2,t) =  mpc.baseMVA*ADMM.var.p_global{t}(2,1); 
    p_inout(3,t) =  mpc.baseMVA*ADMM.var.p_global{t}(3,1); 

    
    q_inout(1,t) =  -ADMM.var.q_global{t}(1,1);
    q_inout(2,t) =  ADMM.var.q_global{t}(2,1); 
    q_inout(3,t) =  ADMM.var.q_global{t}(3,1);   
end

figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;

set(gcf,'Units','inches','Position',figSize);

plot(p_inout(1,:),'-*','DisplayName','MG 1 - bus 7'); hold on;
plot(p_inout(2,:),'-*','DisplayName','MG 2 - bus 7'); hold on;
plot(p_inout(3,:),'-*','DisplayName','MG 3 - bus 7'); hold on;

x = [0, number_of_market_period+2, number_of_market_period+2, 0];
y = [0, 0, 4, 4];
c = [0 0.5 0];
fill(x, y, c, 'FaceAlpha', 0.2,'HandleVisibility','off')
hold on

txt = 'energy import';
text(0.5,3.5,txt,'HorizontalAlignment','left','FontSize',fntsize)
txt = 'energy export';
text(0.5,-3.5,txt,'HorizontalAlignment','left','FontSize',fntsize)
grid on;
set(0,'DefaultAxesFontSize',fntsize);
xlabel('market period');
ylabel('power exchage[MW]');
filename='power';
legend('Location','North','Orientation','vertical');
legend boxoff 
legend show;
ylim([-4 4]);
export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');


%%
figure()
figSize=[2, 2, 6, 1.5];
fntsize=10;
set(0,'DefaultAxesFontSize',fntsize);
set(gcf,'Units','inches','Position',figSize);
plot(ADMM.norm2aux)
set(gca, 'YScale', 'log')
xlabel('iterations');
ylabel('convergence error');
filename='residual';
grid on;
export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');

results_file = ADMM.norm2aux;

mat_file_name = 'rho_1_e_min_3.mat';
save(mat_file_name,'results_file');

%% soc plot 
figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;

plot(0:number_of_market_period,[soc_init(1) soc_res(1,:)],'b-o','LineStyle', '--','DisplayName','central - BESS 1'); hold on;
plot(0:number_of_market_period,[soc_init(2) soc_res(2,:)],'c-o','LineStyle', '--','DisplayName','central - BESS 2'); hold on;
plot(0:number_of_market_period,[soc_init(3) soc_res(3,:)],'k-o','LineStyle', '--','DisplayName','central - BESS 3'); hold on;


plot(0:number_of_market_period,[soc_init(1) ADMM.var.SOC{1}],'b-*','DisplayName','ADMM - BESS 1'); hold on;
plot(0:number_of_market_period,[soc_init(2) ADMM.var.SOC{2}],'c-*','DisplayName','ADMM - BESS 2'); hold on;
plot(0:number_of_market_period,[soc_init(3) ADMM.var.SOC{3}],'k-*','DisplayName','ADMM - BESS 3'); hold on;


legend('Location','NorthEast','Orientation','vertical');
legend boxoff 
legend show;

set(0,'DefaultAxesFontSize',fntsize);
set(gcf,'Units','inches','Position',figSize);
xlabel('market period');
ylabel('SOC');

filename='soc';
export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');

%% converged price 
temp = zeros(1,time_step);
temp2 = zeros(1,time_step);
temp3 = zeros(1,time_step);
for t = 1:time_step
    temp(1,t) = ADMM.var.Lambda{1,t}(1,1);
    temp2(1,t) = ADMM.var.Lambda{2,t}(1,1);
    temp3(1,t) = ADMM.var.Lambda{3,t}(1,1);
end


%%
figure()
figSize=[2, 2, 6, 2.5];
fntsize=10;


%temp_mesh = [temp/mpc.baseMVA;temp2/mpc.baseMVA;temp3/mpc.baseMVA];
%mesh(temp_mesh);
bar(psp_price,'DisplayName','nodal price at PSP'); hold on;
plot(temp,'r-*','DisplayName','MG 1: inport/export at bus 7'); hold on;
plot(temp2,'b-o','DisplayName','MG 2: inport/export at bus 7'); hold on;
plot(temp3,'y-^','DisplayName','MG 3: inport/export at bus 7'); hold on;
%ylim([84 115]);
xlabel('market period');
ylabel('price [S$/MW]');

yyaxis right
plot(load_profile,'b--','DisplayName','normalized load profile');
ylabel('load profile')


set(0,'DefaultAxesFontSize',fntsize);
set(gcf,'Units','inches','Position',figSize);

legend('Location','North','Orientation','vertical');
legend show;
legend boxoff 
filename='import_export_price';
export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');

%%
%%
% figure(1);
% bar(psp_price);
% 
% ylim([80 100]);
% figSize=[2, 2, 6, 1.5];
% fntsize=10;
% set(0,'DefaultAxesFontSize',fntsize);
% set(gcf,'Units','inches','Position',figSize);
% xlabel('market period of 24 market period');
% ylabel('price [S$/MW]');
% grid on;
% hold off;
% 
% 
% filename='psp_price';
% export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');

%%
% figure()
% figSize=[2, 2, 6, 3.5];
% fntsize=10;
% set(0,'DefaultAxesFontSize',fntsize);
% set(gcf,'Units','inches','Position',figSize);
% subplot(1,2,1);
% plot(lambda_rec(1,1:end)/mpc.baseMVA,'y','DisplayName','bus 7, grid 1'); hold on;
% plot(lambda_rec(2,1:end)/mpc.baseMVA,'m','DisplayName','bus 7, grid 2'); hold on;
% plot(lambda_rec(3,1:end)/mpc.baseMVA,'c','DisplayName','bus 7, grid 3'); hold on;
% plot(lambda_rec(4,1:end)/mpc.baseMVA,'r','DisplayName','bus 15, grid 3'); hold on;
% plot(lambda_rec(5,1:end)/mpc.baseMVA,'g','DisplayName','bus 15, grid 4'); hold on;
% plot(lambda_rec(6,1:end)/mpc.baseMVA,'b','DisplayName','bus 24, grid 4'); hold on;
% plot(lambda_rec(7,1:end)/mpc.baseMVA,'k','DisplayName','bus 24, grid 5'); hold on;
% legend('Location','NorthEast','Orientation','vertical');
% legend show;
% xlabel('iterations'); grid on;
% ylabel('active power exchange price');
% hold off;
% 
% subplot(1,2,2);
% plot(kappa_rec(1,1:end)/mpc.baseMVA,'DisplayName','bus 7, grid 1'); hold on;
% plot(kappa_rec(2,1:end)/mpc.baseMVA,'DisplayName','bus 7, grid 2'); hold on;
% plot(kappa_rec(3,1:end)/mpc.baseMVA,'DisplayName','bus 7, grid 3'); hold on;
% plot(kappa_rec(4,1:end)/mpc.baseMVA,'DisplayName','bus 15, grid 3'); hold on;
% plot(kappa_rec(5,1:end)/mpc.baseMVA,'DisplayName','bus 15, grid 4'); hold on;
% plot(kappa_rec(6,1:end)/mpc.baseMVA,'DisplayName','bus 24, grid 4'); hold on;
% plot(kappa_rec(7,1:end)/mpc.baseMVA,'DisplayName','bus 24, grid 5'); hold on;
% legend('Location','NorthEast','Orientation','vertical');
% legend show;
% xlabel('iterations');
% ylabel('reactive power exchange price'); grid on;
% hold off;
% filename='sigma_congested';
% export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');
% 
% 
% %%



% %%
figure()
% plot(obj_v)
% 
% %% voltage 
% figure();
% set(gcf,'Units','inches','Position',figSize);
% set(0,'DefaultAxesFontSize',fntsize)
% set(0, 'DefaultAxesFontName','Times New Roman');
% set(gcf,'color','white')
% 
% acopfres = runopf(mpc);
% 
% plot(acopfres.bus(:,VM)','k-*','DisplayName','ACOPF')
% ylim([0.90 1.1])
% hold on;
%
plot(2:141, v_mag_central_res,'-*'); hold on;
plot(mpc.clusterBus{1}(2:end),ADMM.var.VoltageProfile{1}(:,2)','ko','DisplayName','ADMM-grid 1')
hold on;
plot(mpc.clusterBus{2}(2:end),ADMM.var.VoltageProfile{2}(:,2)','ro','DisplayName','ADMM-grid 2')
hold on;
plot(mpc.clusterBus{3}(2:end),ADMM.var.VoltageProfile{3}(:,2)','bo','DisplayName','ADMM-grid 3')
hold on;
% plot(mpc.clusterBus{4}(2:end),ADMM.var.Voltage{4}','mo','DisplayName','ADMM-grid 4')
% hold on;
% plot(mpc.clusterBus{5}(2:end),ADMM.var.Voltage{5}','go','DisplayName','ADMM-grid 5')
% hold off;
% 
% grid on
% %YMinorGrid on
% %set(gca, 'YScale', 'log')
% %set(gca,'YMinorTick','on')
% box on
% legend('Location','SouthEast','Orientation','horizontal');
% xlabel('Node index');
% ylabel('voltage profile');
% % 
% %% DLMPs 
% 
% mpc_cluster_array(1).gencost(5,5) = 3;
% mpc_cluster_array(1).gen(4,PMAX) = 0.61;
% mpc_cluster_array(1).gen(4,QMIN) = -1.15;
% 
% acopfres2 = runopf(mpc_cluster_array(1));
% 
% 
% DLMP_ACOPF = acopfres2.bus(:,LAM_P);
% DLMP_ACOPFQ = acopfres2.bus(:,LAM_Q);
% 
% figure(3)
% plot(P_DLMPs,'b-*','DisplayName','propsed distributed DLMPs'); hold on;
% plot(DLMP_ACOPF(2:end),'g-o','DisplayName','ACOPF using IPM-Matpower'); hold on;
% ylim([16.5 18.5]);
% figSize=[2, 2, 6, 2.5];
% fntsize=10;
% set(0,'DefaultAxesFontSize',fntsize);
% set(gcf,'Units','inches','Position',figSize);
% xlabel('node index in smart grid 1');
% ylabel('DLMPs for active power in $/MW');
% legend show;
% grid on;
% hold off;
% filename='dlmp_p';
% export_fig(gcf,filename,'-nocrop','-pdf','-r300','-painters');
% Q price does not work well due to the export... reason not clear 
% figure(4) 
% plot(Q_DLMPs,'r-*'); 
% hold on;
% plot(DLMP_ACOPFQ(2:end),'k-o'); 
% ylim([0 5]);
% hold off;
%%
% plot(ADMM.norm2aux)
% set(gca, 'YScale', 'log')
toc