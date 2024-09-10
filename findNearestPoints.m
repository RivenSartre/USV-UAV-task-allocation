function [in, out ]= findNearestPoints(listnoe, weightMatrix)
    % point: 点集
    % list: 抽出点的序号列表
    % weightMatrix: 所有点之间的权重矩阵

    % 获取点集总数
    numPoints = size(weightMatrix, 1);
    
    % 获取不在 list 中的点的索引
    remainingIndices = setdiff(1:numPoints, listnoe);

    % 获取当前点到 list 中所有点的权重距离
    distances = weightMatrix(listnoe,remainingIndices);
    a=min(min(distances));
    % 找到最近的点及其权重
    [in, out]=find(weightMatrix(listnoe,:)==a);
    
    
end