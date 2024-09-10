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