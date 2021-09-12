function [p] = check_parameter_constrains_chem(x)
% check parameter constrains of parameter input vector x, p=1 if fullfilled and 0 if not
% select parameters:
% growth rates higher than dead rates, and dead rate of active higher than
% dormant bacteria (line 1-3);
p=x(2)>=x(9) & x(9)>=x(10)... %bacteria A
    & x(20)>=x(23) & x(23)>=x(24)...%bacteria B
    & x(42)>=x(44) & x(44)>=x(45)...%bacteria D
    & x(59)>max(x([58,60:64]),[],2) & x(58)> max(x([60,62,64]),[],2); %order of sorption coefficients: HY highest; AT>DEA,NE,CY
end