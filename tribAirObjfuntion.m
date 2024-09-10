function [objf] = tribAirObjfuntion(x, AirSet, SurSet, us, CostA)
    % X = [x_1,y_1,x_2,y_2,us,ua];
    X1 = [x(1),x(2)];
    X2 = [x(3),x(4)];
    u = x(5);
    v = x(6);

    Ts1 = norm(SurSet(1,:) - X1)/us;
    Ts2 = norm(X2 - SurSet(2,:))/us;
    Tss = norm(X1 - X2)/u;
    
    Ta1 = norm(AirSet(1,:) - X1)/v;
    Ta2 = norm(X2 - AirSet(3,:))/v;
    Ta12 = norm(AirSet(2,:) - AirSet(1,:))/v;
    Ta23 = norm(AirSet(3,:) - AirSet(2,:))/v;
    Taa1 = CostA(1);
    Taa2 = CostA(2);
    Taa3 = CostA(3);

    lambda1 = 10^(-15);
    lambda2 = 1;
    
    F = Ts1 + Ts2 + Tss;
    Ta = Ta1 + Ta2 + Taa1 + Ta12 + Taa2 + Taa3 + Ta23;
    F1 = lambda1 * Tss * abs(us-u);
    F2 = lambda2 * abs(Tss-Ta);
    objf = F+F1+F2;
end