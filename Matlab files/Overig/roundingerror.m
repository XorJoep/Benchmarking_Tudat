% Graph truncation error
% numericaly calculate the integral of the 3/2*sqrt(x) from x=0 to 1.
clear all;

% exact answer to the question
exactans = 1;
fog=123456789;

%settings of the integrator
resolution = 50;
length = 17*resolution; % do NOT increase the number before the resolution, max 19

result = zeros(1,length);
for i=1:length
    steps = max(exp(i/resolution),2);
    x=linspace(0,1,steps);
    stepsize = x(2);
    result(i,1) = stepsize;
    
    y=(sqrt(x)*(3/2))+fog;
    result(i,2) = (sum(y)-fog*size(y,2))*stepsize;
end
error = abs((exactans-result(:,2)));
loglog(result(:,1),error)
title('Rounding and Truncation error')
xlabel('Stepsize')
ylabel('Relative Error')

saveas(gcf,'ErrorType.png')