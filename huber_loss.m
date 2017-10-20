function [h_loss,huber_all] = huber_loss(raw,processed,sigma)

% subtract the two vectors 
resid = raw - processed;

% build up huber loss for all points
huber_all = (sigma.^2)*(sqrt(1+(resid/sigma).^2)-1);

h_loss = mean(huber_all);

end