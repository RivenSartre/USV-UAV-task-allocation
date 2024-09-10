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