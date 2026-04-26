% This file is a part of the MC2 toolbox developed by Y. Mohammand and T. Nishida.
%Please do not remove this comment
%
% Using this file is governed by the license of MC2 which you can find in LICENSE.md
% 
% You can find more information about this toolbox here:
% - Yasser Mohammad and Toyoaki Nishida, "MC2: An Integrated Toolbox for Change, Causality, 
%   and Motif Discovery", 29th International Conference on Industrial, Engineering & 
%   Other Applications of Applied Intelligent Systems (IEA/AIE) 2016, pp. 128 -- 141.
% - Yasser Mohammad and Toyoaki Nishida, "Data Mining for Social Robotics", Springer 2016.
%

function theta=adjustAngle(theta)
% ensures that the angle is between -pi and pi
% if numel(theta)==1
%    if isinf(theta)
%        return;
%    end
%    if(theta>pi)
%        theta=(theta-floor(theta/pi)*pi)-pi;
%    elseif theta<-pi
%        theta=pi-(-theta-floor(-theta/pi)*pi);
%    end    
% else
%     theta(~isinf(theta) && theta>pi)=(theta(~isinf(theta) && theta>pi)...
%         -floor(theta(~isinf(theta) && theta>pi)/pi)*pi)-pi;
%     
%     theta(~isinf(theta) && theta<-pi)=pi-(-theta(~isinf(theta) && theta<-pi)...
%         -floor(-theta(~isinf(theta) && theta<-pi)/pi)*pi);
    theta=theta+pi;
    theta=theta-(2*pi).*floor(theta./(2*pi));    
    theta=theta-pi;
    theta(isnan(theta))=-inf;    
% end
end