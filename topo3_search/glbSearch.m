clear all
close all

seed = 0
rng default % For reproducibility
rng(seed)

%links and their capacities
load topo35.mat
load('../demandMaxPool.mat') %get demandMaxPool
load('../demandMinPool.mat') %get demandMinPool

numIter = size(demandMaxPool,1);

lambda1 = 0.00133;
beta1 = 0;
lambda2 = 0.08;
beta2 = 350;
lambda3 = 0.03651;
beta3 = 0.5;
lambda4 = 0.00229;
beta4 = 1;

% lambda2_t2 = 0.09
% beta2_t2 = 600

numTraffType = 4;
numSta = size(pathInNode,1);
numPath = size(pathInNode,2);
numFlows = numTraffType * numPath;

fun = @(rate) - sum(sum(rate(1:numFlows/4) * lambda1 + beta1) + ...
                sum(1 ./ (1 + exp(-(rate(numFlows/4 + 1 :numFlows/2) - beta2)* lambda2))) +...
                sum(rate(numFlows/2+1:numFlows*3/4).^ beta3 * lambda3) + ...
                sum(log(rate(numFlows*3/4 + 1 : end)*lambda4 +beta4)));
            
%initialise writing files
demandMaxLog = zeros(numIter,numFlows);
demandMinLog = zeros(numIter,numFlows);
time = zeros(numIter,1);
rateLog = zeros(numIter,numFlows);
utilityTotal = zeros(numIter,1);
utilityIndiv = zeros(numIter,numFlows);
utilityPerType = zeros(numIter,numTraffType);
save('output.mat')

for iter = 1:numIter
    disp('Iteration')
    disp(iter)

    flowDemandMax = demandMaxPool(iter,:);
    flowDemandMin = demandMinPool(iter,:);
    
    %unroll constraints
    time_ub = ones(numSta*2,1);
    pathInOutNode_unroll = [repmat(pathInNode, 1, numTraffType);repmat(pathOutNode, 1, numTraffType)];
    pathInOutCapa_unroll = [repmat(pathInNodeCapacity, 1, numTraffType);repmat(pathOutNodeCapacity, 1, numTraffType)];
    unroll_coef = pathInOutNode_unroll./pathInOutCapa_unroll;

    %bounds of flow rate value
    upper_bound = flowDemandMax  %reshape(flowDemand',numel(flowDemand),1);
    lower_bound = flowDemandMin  %min(upper_bound, 10*randi(10, size(upper_bound)));

    rate0 = flowDemandMin; %zeros(size(flowDemandMax));
                                
    gs = GlobalSearch;

    problem = createOptimProblem('fmincon','x0',rate0,...
    'objective',fun,'lb',lower_bound,'ub',upper_bound,'Aineq',unroll_coef,...
    'bineq',time_ub);

    tic
    [flowRate_opt,utility] = run(gs,problem)
    toc

%   multiSearch
%     rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
%     problem = createOptimProblem('fmincon','objective',...
%     fun,'x0',rate0,'lb',lower_bound,'ub',upper_bound,'Aineq',unroll_coef,...
%     'bineq',time_ub,'options',opts);
% 
%     ms = MultiStart;
%     
%     %optimization
%     tic
%     [flowRate_opt,utility] = run(ms,problem,20)
%     toc

    %flowRate_opt_round = round(flowRate_opt/50)*50;

    utility_indiv= [(flowRate_opt(1:numFlows/4) * lambda1 + beta1), ...
        1 ./ (1 + exp(-(flowRate_opt(numFlows/4 + 1 :numFlows/2) - beta2)* lambda2)),...
        flowRate_opt(numFlows/2+1:numFlows*3/4).^ beta3 * lambda3,...
        log(flowRate_opt(numFlows*3/4 + 1 : end)*lambda4 +beta4)];

    utility_type = sum(reshape(utility_indiv,numPath, numTraffType)',2);
    
   
    %save results 
    demandMaxLog(iter,:) = upper_bound;
    demandMinLog(iter,:) = lower_bound;
    time(iter) = toc;
    rateLog(iter,:) = flowRate_opt;
    utilityTotal(iter) = -utility;
    utilityIndiv(iter,:) = utility_indiv;
    utilityPerType(iter,:) = utility_type;
    
    save('output.mat','demandMaxLog','demandMinLog','time','rateLog', ...
        'utilityTotal','utilityIndiv','utilityPerType','-append')
end
