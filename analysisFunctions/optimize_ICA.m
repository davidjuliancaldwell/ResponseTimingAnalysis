function [x,history,searchdir] = optimize_ICA(data,varargin)

scale_factor = 1000;
stimChans = [];
meanSub = 0;
orderPoly = 2;

% lower bound for scale factor
lb = [1];
ub = [2000];
% starting scale factor 
x0 = [500];


for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'orderPoly'
            orderPoly = varargin{i+1};
        case 'stimChans'
            stimChans = varargin{i+1};
        case 'meanSub'
            meanSub = varargin{i+1};
        case 'examptrial'
            exampTrial = varargin{i+1};
        case 'fs'
            fs = varargin{i+1};
        case 'x0'
            x0 = varargin{i+1};
    end
    
end


history.x = [];
history.fval = [];
searchdir = [];

% using outputfunction example from matlab

figure
% finite difference is for gradient, discrete gradient , + delta_t, looks
% forward and background, this step will affect whether or not gradient is
% good or not, finite difference is either big or small step for finite
% difference "right direction" - pick up noise 

% diffminchange, epsilon, once values dont change , sensitiviy of the
% problem, bigger epsilon, stop sooner , theory may say dont play with time
% step 

% options = optimoptions(@fmincon,'FiniteDifferenceStepSize',5,'Display','iter','PlotFcn',@optimplotx);
%options = optimoptions(@fmincon,'OutputFcn',@outfun,'FiniteDifferenceStepSize',15,'Display','iter');
options = optimoptions(@fmincon,'OutputFcn',@outfun,'diffminchange',5,'Display','iter');

% 
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)ica_train_optimize(x,data,stimChans,fs,meanSub,orderPoly),...
%     x0,[],[],[],[],lb,ub,[],options);

[x] = fmincon(@(x)ica_train_optimize_subtract(x,data,stimChans,fs,meanSub,orderPoly),...
    x0,[],[],[],[],lb,ub,[],options);

function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.gradient'];
           plot(optimValues.iteration,x(1),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
           text(optimValues.iteration+.15,x(1),... 
                num2str(optimValues.iteration));
           title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
 end

% x is the optimal scale factor


end