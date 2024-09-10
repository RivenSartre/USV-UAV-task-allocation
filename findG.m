function [G] = findG(x, AirSet, SurSet, us, CostA, CostS)
   % X = [x_1,y_1,us,ua];
    X = [x(1),x(2)];
    u = x(3);
    v = x(4);

    Ts0 = norm(SurSet(1,:) - X)/us;     % USV到起点的时间
    AirSet = [X; AirSet];
    if  length(AirSet)>1
        for i=1:length(AirSet)-1
            Tapath(i) = norm(AirSet(i,:) - AirSet(i+1,:))/v;
        end
    end
    Tapath=cumsum(Tapath);
    Taa = Tapath(end) + sum(CostA);                  % UAV完成任务的时间

    % 检查时间到Taa的时候，USV在哪里。
    SurSet(1,:)=X;
    for i=1:length(SurSet)-1
        Tspath(2*i-1) = norm(SurSet(i,:) - SurSet(i+1,:))/us;
        if i<length(CostS)
            Tspath(2*i)=CostS(i);
        end
    end
    Tspath=cumsum(Tspath);
    if Taa<Tspath(1)
        % USV还没到第一个点 idx=0;
        idx=0;
    else
        % USV在任务中，奇数是任务中，偶数是路程中。
        idx = find(Tspath < Taa, 1, 'last');
    end
    
    if mod(idx,2)==1 % 在任务中。
        n = 1+ (idx+1)/2; % 在第n个电桩
        Tplus = norm(AirSet(end,:) - SurSet(n,:))/v;
        Tnow = Tplus + Tspath(end);
        G = SurSet(n,:);
    else  % 在路程中
        n = (idx)/2+1 ; % 刚刚离开第n个桩
       % Ttemp=Taa-Tspath(idx);
        
        Snow = SurSet(n,:);
        Anow = AirSet(end,:);

        vec = AirSet(end,:) - SurSet(n,:);
        vec = vec/norm(vec);  % n桩到空点的方向向量
        if isnan(vec(1))
            n = 1+ (idx)/2; % 在第n个电桩
            Tplus = norm(AirSet(end,:) - SurSet(n,:))/v;
            Tnow = Tplus + Tspath(end);
            G = SurSet(n,:);
        else
            %vec = vec/norm(vec);  % n桩到空点的方向向量
            Snow = Snow ;%Ttemp*u *vec;
            k = v/u;
            G = [(Anow(1)+ Snow(1)*k)/(k+1) , (Anow(2)+ Snow(2)*k)/(k+1)];
    
            SurSetplus=[SurSet(1:n,:);G;SurSet(n+1:end,:)];
            for j=1:length(SurSetplus)-1
            Tspath(j) = norm(SurSetplus(j,:) - SurSetplus(j+1,:))/us;
            end
            Tspath=cumsum(Tspath);  
            Tnow = Tspath(end) +sum(CostS);
        end
    end

    lambda1 = 5;
    %lambda2 = 1;
    
    F = Tnow + Ts0;
    F1 = lambda1 * Tnow * abs(us-u);
    %F2 = lambda2 * abs(Tss-Ta);
    objf = F+F1;
end
