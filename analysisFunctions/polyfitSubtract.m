function [dataOut] = polyfitSubtract(dataIn,orderPoly)
% DJC 4-5-2017 - This is a function to subtract off the slowly varying
% component of a signal
%
% input:
%   dataIn: time x channels
%   orderPoly: order of polynomial to use - try 10

dataOut = zeros(size(dataIn));

for i = 1:size(dataIn,2)
    
    data_int = dataIn(:,i);
    [p,s,mu] = polyfit((1:numel(data_int))',data_int,orderPoly);
    f_y = polyval(p,(1:numel(data_int))',[],mu);
    
    dataOut(:,i) = dataIn(:,i) - f_y;
end

end