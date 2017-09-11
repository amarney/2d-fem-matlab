function f = shaperef1dherm(t, iDegree, iDerivative, iShape)
%if iDegree ~= 1
%    disp(['iDegree neq 1'])
%end
if iShape == 1 % match at (-1)
    if iDerivative == 0
        f = (1/4)*(1 - t).^2.*(2 + t);
    elseif iDerivative == 1
        f = (1/4)*(3*t.^2 - 3);
        %f = (1/2)*(-1 + t).*(1 + t) + (1/4)*(1 + t).^2;
    elseif iDerivative == 2
        f = (3/2)*t;
        %f = 1 + (1/2)*(-1 + t) + t;
    elseif iDerivative == 3
        f = (3/2)*ones(size(t));
    elseif iDerivative > 3
        f = zeros(size(t));
    end
elseif iShape == 2 % match at (1)
    if iDerivative == 0
        f = (1/4)*(1 + t).^2.*(2 - t);
    elseif iDerivative == 1
        f = (1/2)*(2 - t).*(1 + t) - (1/4)*(1 + t).^2;
    elseif iDerivative == 2
        f = -1 + (2 - t)/2 - t;
    elseif iDerivative == 3
        f = -(3/2)*ones(size(t));
    elseif iDerivative > 3
        f = zeros(size(t));
    end
elseif iShape == 3 % deriv at (-1)
    if iDerivative == 0
        f = (1/4)*(1 - t).^2.*(t + 1);
    elseif iDerivative == 1
        f = (1/4)*(1 - t).^2 - (1/2)*(1 - t).*(1 + t);
    elseif iDerivative == 2
        f = -1 + t + (1 + t)/2;
    elseif iDerivative == 3
        f = (3/2)*ones(size(t));
    elseif iDerivative > 3
        f = zeros(size(t));
    end
elseif iShape == 4 % deriv at (1)
    if iDerivative == 0
        f = (1/4)*(1 + t).^2.*(t - 1);
    elseif iDerivative == 1
        f = (1/2)*(-1 + t).*(1 + t) + (1/4)*(1 + t).^2;
    elseif iDerivative == 2
        f = 1 + (1/2)*(-1 + t) + t;
    elseif iDerivative == 3
        f = (3/2)*ones(size(t));
    elseif iDerivative > 3
        f = zeros(size(t));
    end
else
    disp(['iShape = ', int2str(iShape)])
    disp('This program is not for this value of iShape')
    disp('Terminate this program by Ctrl+C')
    pause
end
end
