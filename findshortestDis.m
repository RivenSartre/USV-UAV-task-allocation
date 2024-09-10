function short = findshortestDis(node, list, distanceMatrix)
    
    order=distanceMatrix(list,node);

    short=min(order);
    
end