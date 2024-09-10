function [objf] = sinoAirObjfuntionD(x, AirSet, SurSet, us, CostA, CostS)
    % X = [x_1,y_1,x_2,y_2,us,ua];
    X = [x(1),x(2)];
    %G = [x(3),x(4)];
    u = x(3);
    v = x(4);

    Ts0 = norm(SurSet(1,:) - X)/us;     % USV到起点的时间

    Ta1 = norm(AirSet(1,:) - X)/v;
    Taa = CostA + Ta1;                  % UAV完成任务的时间

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
        % USV还没到第一个点
        idx=1;
    elseif Taa>Tspath(end)
        % USV已经到最后的点了
        idx=length(Tspath)+1;
    else
        % USV在任务中，奇数是任务中，偶数是路程中。
        idx = find(Tspath < Taa, 1, 'last');
    end
    
    if mod(idx,2)==1 % 在任务中。
        Tplus = norm(AirSet(end,:) - SurSet(idx+1,:))/v;
        Tnow = Tplus + Tspath(end) + Ts0;
        G = SurSet(idx+1,:);
    else  % 在路程中
        Snow = SurSet((idx)/2+1,:);
        Anow = AirSet(end,:);
        k = v/u;
        G = [(Anow(1)+ Snow(1)*k)*(k+1) , (Anow(2)+ Snow(2)*k)*(k+1)];

        SurSet=[SurSet(1:(idx)/2+1,:);G;SurSet((idx)/2+2:end,:)];
        %CostS = [CostS(1:(idx)/2+1,:);0;CostS((idx)/2+2:end,:)];
        for j=1:length(SurSet)-1
        Tspath(j) = norm(SurSet(j,:) - SurSet(j+1,:))/us;
        end
        Tspath=cumsum(Tspath);  
        Tnow = Ts0 + Tspath(end) +sum(CostS);
    end

    lambda1 = 5;
    %lambda2 = 1;
    
    F = Tnow;
    F1 = lambda1 * Tnow * abs(us-u);
    %F2 = lambda2 * abs(Tss-Ta);
    objf = F+F1;
end
