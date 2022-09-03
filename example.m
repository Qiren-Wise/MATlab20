clear  %多目标优化
clc

MultiObj.fun = @(x) [-10.*(exp(-0.2.*sqrt(x(:,1).^2 + x(:,2).^2)) + exp(-0.2.*sqrt(x(:,2).^2 + x(:,3)))...
               sum(abs(x).^0.8 + 5.*sin(x.^3),2)];
MultiObj.nVar = 3; % 3目标
MultiObj.var_min = -5.*ones(1,MultiObj.nVar); %目标范围  位置
MultiObj.var_max = 5.*ones(1,MultiObj.nVar);

params.Np = 200; %种群个数
params.Nr = 150;  %前沿放置点数
params.maxgen = 100;  %迭代次数
params.W = 0.4; % W C1 C2 三个参数
params.C1 = 2;
params.C2 = 2;
params.ngrid = 20; %网格个数
params.maxvel = 5; %最大速度
params.u_mut = 0.5; %变异率

REP = MOPSO(params,MultiObj); 
