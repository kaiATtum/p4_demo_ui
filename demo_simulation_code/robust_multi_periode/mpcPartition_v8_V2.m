function [ mpc_cluster_array ] = mpcPartition_v8_V2( mpc )
%UNTITLED Summary of this function goes here
% v8: include the energy storage
%   create sub mpcs under current mpc according to clustering formation 
    define_constants;
    % define struct array type 
    % default 
    mpc_cluster_array.version = '2';
    mpc_cluster_array.baseMVA = 10;
    % ajustable 
    mpc_cluster_array.bus = [];
    mpc_cluster_array.gen = [];
    mpc_cluster_array.branch = [];
    mpc_cluster_array.gencost = [];
    mpc_cluster_array.VA_slack = [];
    mpc_cluster_array.storage = [];
    
    for c = 1:mpc.cluster.number
        bus_temp = mpc.bus(mpc.clusterBus{c},:);
        bus_temp(:,1) = (1:mpc.cluster.nbu{c})';
        
        mpc_cluster_array(c).bus = bus_temp;        
        mpc_cluster_array(c).gen = mpc.clusterGen{c};
        mpc_cluster_array(c).gencost = mpc.clusterGenCost{c};
        mpc_cluster_array(c).storage = mpc.clusterStorage{c};
        mpc_cluster_array(c).branch = mpc.clusterBranch{c};
        % default
        mpc_cluster_array(c).version = '2';
        mpc_cluster_array(c).baseMVA = mpc_cluster_array.baseMVA;
        mpc_cluster_array(c).VA_slack = 0;
             
        
        if c == 1 % add coupling bus as dg bus for cluster 1 
            mpc_cluster_array(c).gen = [mpc_cluster_array(c).gen; mpc_cluster_array(c).gen(1,:)];
            mpc_cluster_array(c).gen(4,1) = 7; % define the coupled bus index 
            
            %add price for coupled bus (become the dg now)
            k = size(mpc_cluster_array(c).gen,1)-1;%number of DGs 
            mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost(1:k,:); mpc_cluster_array(c).gencost(1,:); mpc_cluster_array(c).gencost(k+1:end,:)];% for P  % value doesnot matter           
            mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost; mpc_cluster_array(c).gencost(k+2,:)]; %for Q % value doesnot matter
            
            % relax the limit on coupled bus 
            mpc_cluster_array(c).gen(4,PMIN) = - 200;
            mpc_cluster_array(c).gen(4,PMAX) = 200;
            mpc_cluster_array(c).gen(4,QMIN) = - 200;
            mpc_cluster_array(c).gen(4,QMAX) = 200;
            
        else
            % add slack bus as additional gen for children bus 
            rootIdx = mpc.clusterBus{c}(1);
            rootIdxIternal = find(mpc.clusterBus{c}==rootIdx);
            % add one more gen
            mpc_cluster_array(c).gen = [mpc_cluster_array(c).gen(1,:); mpc_cluster_array(c).gen];
            mpc_cluster_array(c).gen(1,1) = rootIdxIternal;
            
            % change bus type 
            mpc_cluster_array(c).bus(rootIdxIternal,BUS_TYPE) = 3;% slack bus 
            
            % add one more gencost,  gencost value does not matter (only because pf res needed)
            % add one more price for P and Q respectively 
            k = size(mpc_cluster_array(c).gen,1)-1;%number of DGs 
            mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost(1:k,:); mpc_cluster_array(c).gencost(k+1,:); mpc_cluster_array(c).gencost(k+1:end,:)];% for Q
            mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost(1,:); mpc_cluster_array(c).gencost]; %for P    
            
            % relax the limit on coupled bus            
            mpc_cluster_array(c).gen(1,PMIN) = - 200;
            mpc_cluster_array(c).gen(1,PMAX) = 200;
            mpc_cluster_array(c).gen(1,QMIN) = - 200;
            mpc_cluster_array(c).gen(1,QMAX) = 200;
                       
            % for region 3-4 add one more bus for coupling 
%             if c == 3||c == 4
%                 % find children coupled bus index 
%                 if c == 3
%                 rootIdx = 15;
%                 elseif c == 4
%                 rootIdx = 24;
%                 end
%                 rootIdxIternal = find(mpc.clusterBus{c}==rootIdx);    
% 
%                 mpc_cluster_array(c).gen = [mpc_cluster_array(c).gen; mpc_cluster_array(c).gen(1,:)];
%                 mpc_cluster_array(c).gen(3,1) = rootIdxIternal; % define the coupled bus index 
% 
%                 %add price for coupled bus (become the dg now)
%                 k = size(mpc_cluster_array(c).gen,1)-1;%number of DGs 
%                 mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost(1:k,:); mpc_cluster_array(c).gencost(1,:); mpc_cluster_array(c).gencost(k+1:end,:)];% for P  % value doesnot matter           
%                 mpc_cluster_array(c).gencost = [mpc_cluster_array(c).gencost; mpc_cluster_array(c).gencost(k+2,:)]; %for Q % value doesnot matter
% 
%                 % relax the limit on coupled bus 
%                 mpc_cluster_array(c).gen(3,PMIN) = - 200;
%                 mpc_cluster_array(c).gen(3,PMAX) = 200;
%                 mpc_cluster_array(c).gen(3,QMIN) = - 200;
%                 mpc_cluster_array(c).gen(3,QMAX) = 200;
%                 
%                 if c == 3
%                 	mpc_cluster_array(c).gen(1,PMIN) = - 0.5;
%                     mpc_cluster_array(c).gen(1,PMAX) = 0.5;
%                     mpc_cluster_array(c).gen(1,QMIN) = - 2;
%                     mpc_cluster_array(c).gen(1,QMAX) = -2;
%                 end
%             end
        
        end

end

