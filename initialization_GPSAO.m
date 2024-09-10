
function pop = initialization_GPSAO(SearchAgents_no,dim,ub,lb)
%
p = zeros(SearchAgents_no,dim);
prime_number_min = dim*2 +3;
% 
while 1
    if isprime(prime_number_min)==1
        break;
    else
       prime_number_min = prime_number_min + 1;
    end
end

for i = 1:SearchAgents_no
    for j = 1:dim
        r = mod(2*cos(2*pi*j/prime_number_min)*i,1);% 
        p(i,j) = lb(j)+r*(ub(j)-lb(j));
    end
end
pop = p;
end