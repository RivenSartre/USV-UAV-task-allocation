%function [Shortest_Length, Shortest_Route, Gbest_fitness, Length_ave] = PSOTSPsolve(citys, genmax, m)

    %% 1.读入数据
    citys = rand(10, 2)*40; % 随机生成点作为示例
    
    %% 2.计算距离矩阵
    n=size(citys,1);
    D=zeros(n,n);
    for i=1:n
        for j=i+1:n
            D(i,j)=sqrt(sum((citys(i,:)-citys(j,:)).^2));
            D(j,i)=D(i,j);
        end
    end
    
    %% 3.初始化参数
    c1=0.1;                         % 个体学习因子
    c2=0.075;                       % 社会学习因子
    w=1;                            % 惯性因子
    m=200;                          % 粒子数量
    pop=zeros(m,n);                 % 粒子位置
    v=zeros(m,n);                   % 粒子速度
    gen=1;                          % 迭代计数器
    genmax=300;                     % 迭代次数
    fitness=zeros(m,1);             % 适应度函数值
    Pbest=zeros(m,n);               % 个体极值路径
    Pbest_fitness=zeros(m,1);       % 个体极值
    Gbest=zeros(genmax,n);          % 群体极值路径
    Gbest_fitness=zeros(genmax,1);  % 群体极值
    Length_ave=zeros(genmax,1);     % 各代路径的平均长度
    ws=1;                           % 惯性因子最大值
    we=0.8;                         % 惯性因子最小值
    
    %% 4.产生初始粒子
    % 4.1随机产生粒子初始位置和速度
    for i=1:m
        pop(i,:)=randperm(n);
        v(i,:)=randperm(n);
    end
    
    % 4.2计算粒子适应度函数值
    for i=1:m
        for j=1:n-1
            fitness(i)=fitness(i) + D(pop(i,j),pop(i,j+1));
        end
        fitness(i)=fitness(i) + D(pop(i,end),pop(i,1));
    
        %fitness(i) = PNGfit(citys,pop(i,:));
    end
    
    % 4.3计算个体极值和群体极值
    Pbest_fitness=fitness;
    Pbest=pop;
    [Gbest_fitness(1),min_index]=min(fitness);
    Gbest(1,:)=pop(min_index,:);
    Length_ave(1)=mean(fitness);
    
    tic
    % 添加进度条
    h = waitbar(0, ['已完成:0%   搜索中...   已用时:', num2str(toc)]);
    WaitbarInter = genmax / 100;  
    %% 5.迭代寻优
    while gen<genmax
        % 5.1更新迭代次数与惯性因子
        gen=gen+1;
        w = ws - (ws-we)*(gen/genmax)^2;
    
        % 5.2更新速度
        %个体极值修正部分
        change1=position_minus_position(Pbest,pop);
        change1=constant_times_velocity(c1,change1);
        %群体极值修正部分
        change2=position_minus_position(repmat(Gbest(gen-1,:),m,1),pop);
        change2=constant_times_velocity(c2,change2);
        %原速度部分
        v=constant_times_velocity(w,v);
        %修正速度
        for i=1:m
            for j=1:n
                if change1(i,j)~=0
                    v(i,j)=change1(i,j);
                end
                if change2(i,j)~=0
                    v(i,j)=change2(i,j);
                end
            end
        end
    
        % 5.3更新位置
        pop=position_plus_velocity(pop,v);
    
        % 5.4适应度函数值更新
        fitness=zeros(m,1);
        for i=1:m
            for j=1:n-1
                fitness(i)=fitness(i) + D(pop(i,j),pop(i,j+1));
            end
            fitness(i)=fitness(i) + D(pop(i,end),pop(i,1));
            %fitness(i) = PNGfit(citys,pop(i,:));
        end
    
        % 5.5个体极值与群体极值更新
        for i=1:m
            if fitness(i)<Pbest_fitness(i)
                Pbest_fitness(i)=fitness(i);
                Pbest(i,:)=pop(i,:);
            end
        end
    
        [minvalue,min_index]=min(fitness);
        if minvalue<Gbest_fitness(gen-1)
            Gbest_fitness(gen)=minvalue;
            Gbest(gen,:)=pop(min_index,:);
        else
            Gbest_fitness(gen)=Gbest_fitness(gen-1);
            Gbest(gen,:)=Gbest(gen-1,:);
        end
        Length_ave(gen)=mean(fitness);
    
        % 生成进度条
        if mod(gen, WaitbarInter) == 0
            waitbar(gen / genmax, h, ['已完成:', num2str(gen / genmax * 100), '%   搜索中...   已用时:', num2str(toc)]);
        end
    end

    [Shortest_Length,index] = min(Gbest_fitness);
    Shortest_Route = Gbest(index,:);

%end





%% 辅助函数
function change=position_minus_position(best,pop)
%记录将pop变成best的交换序列
for i=1:size(best,1)
    for j=1:size(best,2)
        change(i,j)=find(pop(i,:)==best(i,j));
        temp=pop(i,j);
        pop(i,j)=pop(i,change(i,j));
        pop(i,change(i,j))=temp;
    end
end
end
function change = constant_times_velocity(constant,change)
% 以一定概率保留交换序列
for i=1:size(change,1)
    for j=1:size(change,2)
        if rand>constant
            change(i,j)=0;
        end
    end
end
end
function pop = position_plus_velocity(pop,v)
%利用速度记录的交换序列进行位置修正
for i=1:size(pop,1)
    for j=1:size(pop,2)
        if v(i,j)~=0
            temp=pop(i,j);
            pop(i,j)=pop(i,v(i,j));
            pop(i,v(i,j))=temp;
        end
    end
end
end