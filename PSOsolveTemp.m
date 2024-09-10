function [GBest, GBestFitness ,GBestLand] = PSOsolveTemp(Model2, Vel, ub, lb, iter_Max, N)
D = 2; % 决策变量维度

V_Max = (ub - lb) * 0.15; 
V_Min = -V_Max; % 速度上下限
w_Max = 0.9; w_Min = 0.4; % 惯性因子上下限
c1 = 1.5; c2 = 1.5; % 两个学习因子

WaitbarInter = iter_Max / 100; % 一个和进度条有关的参数

% 粒子初始化
X = lb + (ub - lb) .* rand(N, D);
V = V_Min + (V_Max - V_Min) .* rand(N, D);

% 写入当前全局最优解和个体最优解
GBestFitness = zeros(1, iter_Max + 1);
PBestFitness = zeros(N, 1);
PBest = X;
PBestLand = zeros(N,2);
for i = 1:N
    [PBestFitness(i),PBestLand(i,:)] = Fitness(Model2, Vel, X(i, :));
end

[GBestFitness(1), GBestIndex] = min(PBestFitness);
GBest = PBest(GBestIndex, :);

tic
h = waitbar(0, ['已完成:0%   运算中...用时:', num2str(toc)]);

for iter = 1:iter_Max
    % 初始化当前迭代次数下的最优解
    CurrentIterGbest = GBest;
    CurrentIterGbestFitness = GBestFitness(iter);

    k = 0.6; % 控制因子，一般取0.6
    w = (w_Max - w_Min) * tan(0.875 * (1 - (iter / iter_Max)^k)) + w_Min;

    for i = 1:N
        % 更新个体的位置和速度
        V(i, :) = w * V(i, :) + c1 * rand * (PBest(i, :) - X(i, :)) + c2 * rand * (CurrentIterGbest - X(i, :));
        X(i, :) = X(i, :) + V(i, :);
        % 边界条件限制
        [X(i, :), V(i, :)] = BoundaryLimit(X(i, :), V(i, :), ub, lb, V_Max, V_Min);
        % 计算当前迭代次数下当前粒子的适应度
        [CurrentFitness, CurrentLand]= Fitness(Model2, Vel, X(i, :));
        % 更新当前迭代次数下该粒子的最优解！！！
        % 事实证明，粒子个体最优解的更新应该在一个粒子位置苏都更新完后立即更新！！！
        % 若在当前迭代次数下所有粒子的位置和速度更新完毕后，再更新个体最优解，那么及其难以收敛！！！
        if CurrentFitness < PBestFitness(i)
            PBestFitness(i) = CurrentFitness;
            PBestLand(i,:) = CurrentLand;
            PBest(i, :) = X(i, :);
            % 更新当前迭代次数下的最优解
            if CurrentIterGbestFitness > CurrentFitness
                CurrentIterGbestFitness = CurrentFitness;
                CurrentIterGbestLand = CurrentLand;
                CurrentIterGbest = X(i, :);
            end

        end

    end

    % 依据当前迭代次数下的最优解更新全局最优解
    GBest = CurrentIterGbest;
    GBestFitness(iter + 1) = CurrentIterGbestFitness;
    GBestLand = CurrentIterGbestLand;
    % 展示进度条
    if mod(iter, WaitbarInter) == 0
        waitbar(iter / iter_Max, h, ['已完成:' num2str(iter / iter_Max * 100) ...
        '%   运算中...用时:', num2str(toc),'/',num2str(toc/(iter / iter_Max))])
    end

end
%[Temp_L,~,~ ]= CalcuTempLand(Model1,Vel,GBest);
close(h)

end

% 目标函数
function [CostS, Temp_L]= Fitness(Model2, Vel, X)
    Temp_T = X;
    [Temp_L,place,Flag]= CalcuTempLand(Model2,Vel,Temp_T);

    if Flag==1
        if place==1
            Surnode = [Model2.Surnode(1,:); Temp_T; Temp_L; Model2.Surnode(2:end,:)];
        else 
            Surnode = [Model2.Surnode(1,:); Temp_T; Model2.Surnode(2:place-1,:); Temp_L; Model2.Surnode(place:end,:)];
        end
    elseif Flag==2
        Surnode = [Model2.Surnode(1,:); Temp_T; Model2.Surnode(2:end,:)];
    elseif Flag==0
        Surnode = [Model2.Surnode(1,:); Temp_T; Model2.Surnode(2:end,:); Temp_L;];
    end
    %Surnode = [Model.Surnode(1,:); Temp_T; Temp_L; Model.Surnode(2,:)];
    Airnode = [Temp_T; Model2.Airnode; Temp_L];
    SurLength = 0;
    AirLength = 0;
    us = Vel.us;
    ua = Vel.ua;
    for i=2:length(Surnode)
        SurLength = SurLength + sqrt((Surnode(i,1)-Surnode(i-1,1))^2 + (Surnode(i,2)-Surnode(i-1,2))^2) ;
    end
    for i=2:length(Airnode)
        AirLength = AirLength + sqrt((Airnode(i,1)-Airnode(i-1,1))^2 + (Airnode(i,2)-Airnode(i-1,2))^2) ;
    end

    CostS = SurLength /us;
    CostT = SurLength /us;

end

function [Temp_L,place,Flag ]= CalcuTempLand(Model2,Vel,Temp_T)
    ua = Vel.ua;
    Airnode = [Temp_T; Model2.Airnode];
    AirLength = 0;
    for i=2:length(Airnode)
        AirLength = AirLength + sqrt((Airnode(i,1)-Airnode(i-1,1))^2 + (Airnode(i,2)-Airnode(i-1,2))^2);
    end
    TimeA = sum(Model2.CostAir) + AirLength/ua;


    us = Vel.us;
    Surnode = [Temp_T; Model2.Surnode(2:end,:)];
    CostSur = [0; Model2.CostSur(2:end,:)];
    TimeS(1) = 0;
    item=2;
    Temp_L = [];
    for i=2:length(Surnode)
        TimeS(item) = TimeS(item-1) + sqrt((Surnode(i,1)-Surnode(i-1,1))^2 + (Surnode(i,2)-Surnode(i-1,2))^2) / us;
        if TimeS(item) > TimeA
            % 在航行中UAV完工了
            DeltaV = sqrt((Surnode(i-1,1)-Airnode(end,1))^2 + (Surnode(i-1,2)-Airnode(end,2))^2);
            DeltaT = (DeltaV-(TimeA-TimeS(item-1))*us)/(ua+us);
            dal = DeltaT * ua;
            Temp_L = dal/DeltaV * Airnode(end,:) + (1-dal/DeltaV) * Surnode(i-1,:);

            place = i-1; % 插在Temp_T后i个
            Flag = 1;
            break
        end
        item = item + 1;
        TimeS(item) = TimeS(item-1) + CostSur(i);
        if TimeS(item) > TimeA
            % 在工作中UAV完工了
            Temp_L = Surnode(i,:);

            place = i; % 不用插
            Flag = 2;
            break
        end
        item = item + 1;
    end
    if isempty(Temp_L)
        % 要去最后一个点咯
        DeltaV = sqrt((Surnode(i-1,1)-Airnode(end,1))^2 + (Surnode(i-1,2)-Airnode(end,2))^2);
        DeltaT = (DeltaV-(TimeA-TimeS(item-1))*us)/(ua+us);
        dal = DeltaT * ua;
        Temp_L = dal/DeltaV * Airnode(end,:) + (1-dal/DeltaV) * Surnode(i-1,:);
        %Temp = Surnode(end,:);
        place = length(Surnode);
        Flag = 0;
    end
end

% 边界限制函数
function [result_X, result_V] = BoundaryLimit(X, V, X_Max, X_Min, V_Max, V_Min)

    for i_temp = 1:length(X( 1,:))
        
        if X(i_temp) > X_Max(i_temp)
            X(i_temp) = X_Max(i_temp);
        end

        if X(i_temp) < X_Min(i_temp)
            X(i_temp) = X_Min(i_temp);
        end

    end

    for i_temp = 1:length(V( 1,:))

        if V(i_temp) > V_Max(i_temp)
            V(i_temp) = V_Max(i_temp);
        end

        if V(i_temp) < V_Min(i_temp)
            V(i_temp) = V_Min(i_temp);
        end

    end

    result_X = X;
    result_V = V;

end