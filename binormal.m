function Number = binormal(N,p)
%
% Evaluacion rápida de la binormal 
%
li = length(N);
Number = zeros(size(N));

for j=1:li
if N(j) == 0
    Number(j) = 0;
elseif N(j)<1000
    Number(j) = sum(rand(N(j),1) < p(j));
else
    Number(j) = floor(normrnd(N(j)*p(j),sqrt(N(j)*p(j)*(1-p(j)))));
end
end;
