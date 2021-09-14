function [BLsubtracted] = subtructBL(curve,polyFitOrder,timeResolution)
%subtract growth component from nucleoid length curve

timePoints = transpose(0:timeResolution:(timeResolution*size(curve,1)-1)); % in sec
[p,~,mu] = polyfit(timePoints,curve,polyFitOrder);
baseline = polyval(p,timePoints,[],mu);
BLsubtracted = curve - baseline;
end

