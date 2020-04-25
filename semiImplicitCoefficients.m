
numberOfYEvals = 4;
numberOfFIEvals = 1;
numberOfFEEvals = 3;
numberOfCoeff = numberOfYEvals + numberOfFIEvals + numberOfFEEvals;

Order = 3;

M = zeros(2*Order+1,numberOfCoeff);
b = zeros(2*Order+1,1);

M(1,1:numberOfYEvals) = 1;

for i = 0:numberOfYEvals-1
    for j = 1:Order
        M(j+1,i+1) = (i - 1)^j / factorial(j);
    end
end

for i = 0:numberOfFIEvals-1
    for j = 1:Order
        M(j+1,numberOfYEvals+i+1) = (i - 1)^(j-1) / factorial(j-1);
    end
end

for i = 0:numberOfYEvals-1
    for j = 1:Order
        M(Order + j+1,i+1) = (i - 1)^j / factorial(j);
    end
end

for i = 1:numberOfFEEvals
    for j = 1:Order
        M(Order + j+1,numberOfYEvals+numberOfFIEvals+i) = (i - 1)^(j-1) / factorial(j-1);
    end
end

numberOfExtraConditions = 1;

N = zeros(numberOfExtraConditions,numberOfCoeff);
c = zeros(numberOfExtraConditions,1);

N(1,1) = 1;
c(1) = 1;

% theta = 1;
% N(2,1) = theta;
% N(2,numberOfYEvals+1) = -1;



A = [M; N];
d = [b; c];

size(A)
A
d

coeff = A\d

yCoeff = coeff(1:numberOfYEvals)
implicitCoeff = coeff(numberOfYEvals+1:numberOfYEvals+numberOfFIEvals)
explicitCoeff = coeff(numberOfYEvals+numberOfFIEvals+1:numberOfCoeff)
