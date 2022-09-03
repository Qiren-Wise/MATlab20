function REP = MOPSO(params, MultiObj)   %多目标优化
Np = params.Np;
Nr = params.Nr;
maxgen = params.maxgen;
W = params.W;
C1 = params.C1;
C2 = params.C2;
ngrid = params.ngrid;
maxvel = params.maxvel;
u_mut = params.u_mut;
fun = MultiObj.fun;  %优化
nVar = MultiObj.nVar;
var_min = MultiObj.var_min(:);
var_max = MultiObj.var_max(:);

POS = repmat((var_max-var_min)',Np,1).*rand(Np,nVar)+repamt(var_min',Np,1);
VEL = zeros(Np,nVar);
POS_fit = fun(POS); %适应度情况
PBEST = POS;
PBEST_fit = POS_fit;
DOMINATED = checkDomination(POS_fit); %支配
REP.pos = POS(~DOMINATED,:); %支配个体赋值
REP.pos_fit = POS_fit(~DOMINATED,:); %支配适应度函数赋值
REP = updateGrid(REP,ngrid);    %种群更新网格
maxvel = (var_max - var_min).*maxvel./100;
gen = 1;
display(['Generation #0 - Repossitry size:' num2str(size(ERP.pos,1))]);

stopCondition = false;
while ~stopCondition
    h = selectLeader(REP);
    VEL = W.*VEL + C1 * rand(Np,nVar) .* (PBEST -  POS)...
        +C2 * rand(Np,nVar).*(repmat(REP.pos(h,:),Np,1) - POS);
    POS = POS + VEL;
    POS = mutation(POS,gen,maxgen,Np,var_max,var_min,nVar,u_mut);
    [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min);
    POS_fit = fun(POS);
    
    REP = POS + VEL;
    POS = mutation(POS,gen,maxgen,Np,var_max,var_min,nVar,u_mut);  %变异
    [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min);
    POS_fit = fun(POS);
    
    REP = updateRepository(REP,POS,POS_fit,ngrid); %更新记忆
    if (size(REP.pos,1) > Nr)
        REP = deleteFromRepository(REP,size(REP.pos,1)-Nr,ngrid);
    end
    pos_best = dominates(POS_fit,PBEST_fit);  %检测支配
    best_pos = dominates(PBEST_fit,POS_fit);
    best_pos(rand(Np,1) >=0.5) = 0; %一定概率
    if (sum(pos_best)>1)
        PBEST_fit(pos_best,:)  = POS_fit(pos_best,:); %替换
        PBEST(pos_best,:) = POS(pos_best,:);
        
    end
if (sum(best_pos) >1)
    PBEST_fit(best_pos,:) = POS_fit(best_pos,:);
    PBEST(best_pos,:) = POS(best_pos,:);
end

display(['Generation #' num2str(gen) ' - Repository size:' num2str(size(REP.pos,1))]);

gen = gen + 1;
if (gen > maxgen)
   stopCondition = true;
end
end
     h_rep = plot(REP.posfit(:,1),REP.pos_fit(:,2),'ok');
end

function REP = updateRepository(REP, POS, POS_fit, ngrid)
    DOMINATED = checkDomination(POS_fit);
    REP.pos = [REP.pos;POS(~DOMINATED, :)];
    REP.pos_fit = [REP.pos_fit;POS_fit(~DOMINATED, :)]; %当前状况
    DOMINATED = checkDomination(REP.pos_fit);
    REP.pos_fit = REP.pos_fit(~DOMINATED,:);
    REP.pos = REP.pos(~DOMINATED,:);
    REP = updateGrid(REP,ngrid);
end



function dom_vector = checkDomination(fitness)
Np = size(fitness,1);
dom_vector = zeros(Np,1);%全零值
all_perm = nchoosek(1:Np,2);
all_perm = [all_perm;[all_perm(:,2),all_perm(:,1)]]; %一二列互换

d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
dominated_particles = unique(all_perm(d == 1,2));
dom_vector(dominated_particles) = 1;
end

function d = dominates(x,y)
        d = all(x <=y,2) & any(x < y,2);  %满足 x 支配 y
end

function REP = deleteFromRepository(REP,n_extra,ngrid)
   crowding = zeros(size(REP.pos,1),1);
   for m = 1 : 1 : size(REP.pos_fit,2)
      [m_fit,idx] = sort(REP.pos_fit(:,m),'ascend'); %目标函数值 升序排列
      m_up = [m_fit(2 : end);Inf];
      distance = (m_up - m_down)./(max(m_fit)-min(m_fit));
      [~,idx] = sort(idx,'ascend');
      crowding = crowding + distance(idx);
   end
   crowding(isnan(crowding)) = Inf;
   
   [~,del_idx] = sort(crowding,'ascend');
   del_idx = del_idx(1 : n_extra);  %距离最近的几个数被删除
   REP.pos(del_idx,:) = [];
   REP.pos_fit(del_idx,:) = [];
   REP = updateGrid(REP,ngrid);
end

function selected = selectLeader(REP)
prob = cumsum(REP.quality(:,2));
sel_hyp = REP.quality(find(rand(1,1)*max(prob) <=prob,1,'first'),1);  %选择网格  网格里数越多值越小，概率低
idx = 1 : 1 : length(REP.grid_idx);
selected = idx(REP.grid_idx==sel_hyp);
selected = selected(randi(length(selected)));
  
end

    function REP = updateGrid(REP,ngrid)
         ndim = size(REP.pos_fit,2); %2维
         REP.hypercube_limits = zeros(nrid + 1,ndim);
         for dim = 1:ndim
             REP.hypercube_limits(:,dim) = linspace(min(REP.pos_fit(:,dim)),max(REP.pos_fit(:,dim)),ngrid + 1)';
         end
         
         npar = size(REP.pos_fit,1);
         REP.grid_idx = zeros(npar,1);
         REP.grid_subidx = zeros(npar,ndim);
         for n = 1 : 1 : npar
             idnames = [];
             for d = 1 : 1 : ndim
                 REP.grid_subidx(n,d) = find(REP.pos_fit(n,d) <= REP.hypercube_limits(:,d)',1,'first') - 1; %第一次小于某条线
                 if  (REP.grid_subidx(n,d) ==0)
                     REP.grid_subidx(n,d) = 1;
                 end
                 idnames = [idnames ',' num2str(REP.grid_subidx(n,d))];
             end
             REP.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames'):']);
         end
         REP.quality = zeros(ngrid,2);
         ids = unique(REP.grid_idx);
         for i = 1:length(ids)
             REP.quality(i,1) = ids(i);
             REP.quality(i,2) = 10/sum(REP.gris_idx == ids(i)); %求得网格质量
         end
    end
    
    function [POS,VEL] = checkBoundaries(POS,VEL, maxvel, var_max, var_min)
        Np = size(POS,1);
        MAXLIM = repmat(var_max(:)', Np,1);
        MINLIM = repmat(var_min(:)', Np,1);
        MAXVEL = repmat(maxvel(:)', Np,1);
        MINVEL = repmat(-maxvel(:)', Np,1);
          
           
           VEL(VEL > MAXVEL) = MAXVEL(VEL > MAXVEL);
           VEL(VEL < MINVEL) = MINVEL(VEL < MINVEL);
           VEL(POS > MAXLIM) = (-1).*VEL(POS > MAXLIM);
           POS(POS > MAXLIM) = MAXLIM(POS > MAXLIM);
           VEL(POS < MINLIM) =(-1).*VEL(POS < MINLIM);
           POS(POS < MINLIM) = MINLIM(POS < MINLIM);
    end