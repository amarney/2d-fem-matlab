function f = shape_function_generator_1d_lagrange_reference(t, ...
    iDegree, iDerivative, iShape)
% %--------------------------------------------------------------------------
% %shape_function_generator_1d_lagrange_reference:
% %   returns the values (or derivative values) for the chosen local
% %   finite-element shape function over the reference element [-1,1]. The
% %   shape functions themselves are Lagrange polynomials. Designed to be
% %   called within the function ``shape_function_generator_1d_lagrange''. 
% %
% %   Parameters  :   t - an array from -1 to 1 consisting of coordinates
% %                       where the shape function is being evaluated.
% %                   iShape - index specifying which local
% %                             shape function is evaluated.
% %                   iDerivative - index specifying which
% %                                  derivative of the shape function 
% %                                  is evaluated.
% %                   iDegree - index specifying the degree of the
% %                            Lagrange polynomial being evaluated.
% %
% %   Return      :   f - an array consisting of values (or derivative
% %                       values) corresponding to the coordinates of t.
% %   
% %   Example     :
% %       To plot the third shape function on [-1,1] for second degree
% %       polynomials do
% %
% %       t = -1:0.01:1; iDegree = 2; iDerivative = 0; iShape = 3; 
% %       f = shape_function_generator_1d_lagrange_reference(t, ... 
% %       iDegree, iDerivative, iShape);
% %       plot(t,f)
% %
% %       To plot all the shape functions on [-1,1] for third degree
% %       polynomials do
% %
% %       t = -1:0.01:1; iDegree = 3; iDerivative = 0;
% %       f1 = shape_function_generator_1d_lagrange_reference(t, ... 
% %       iDegree, iDerivative, 1);
% %       f2 = shape_function_generator_1d_lagrange_reference(t, ... 
% %       iDegree, iDerivative, 2);
% %       f3 = shape_function_generator_1d_lagrange_reference(t, ... 
% %       iDegree, iDerivative, 3);
% %       f4 = shape_function_generator_1d_lagrange_reference(t, ... 
% %       iDegree, iDerivative, 4);
% %       plot(t,f1,t,f2,t,f3,t,f4)
% %
% %       To plot the second derivatives, just do
% %       iDerivative = 2;
% %
% %   NOTE: Only goes up to third degree!
% %--------------------------------------------------------------------------
if iDegree == 0
    if iShape == 1
        if iDerivative == 0
            f = ones(size(t));
        elseif iDerivative > 0
            f = zeros(size(t));
        end
    else
        disp(['iShape = ', int2str(iShape)])
        disp('This program is not for this value of iShape')
        disp('Terminate this program by Ctrl+C')
        pause
    end
elseif iDegree == 1
    if iShape == 1
        if iDerivative == 0
            f = (1-t)/2;
        elseif iDerivative == 1
            f = (-1/2)*ones(size(t));
        elseif iDerivative > 1
            f = zeros(size(t));
        end
    elseif iShape == 2
        if iDerivative == 0
           f = (1+t)/2;
        elseif iDerivative == 1
           f = (1/2)*ones(size(t));
        elseif iDerivative > 1
            f = zeros(size(t));
        end
    else
        disp(['iShape = ', int2str(iShape)])
        disp('This program is not for this value of iShape')
        disp('Terminate this program by Ctrl+C')
        pause
    end   
elseif iDegree == 2
    if iShape == 1
        if iDerivative == 0
            f = -(1/2)*(1-t).*t;
        elseif iDerivative == 1
            f = -1/2+t;
        elseif iDerivative == 2
            f = ones(size(t));
        elseif iDerivative > 2
            f = zeros(size(t));
        end
    elseif iShape == 2
        if iDerivative == 0
            f = 1-t.^2;
        elseif iDerivative == 1
            f = -2*t;
        elseif iDerivative == 2
            f = -2*ones(size(t));
        elseif iDerivative > 2
            f = zeros(size(t));
        end
    elseif iShape == 3
        if iDerivative == 0
            f = (1/2)*t.*(1+t);
        elseif iDerivative == 1
            f = 1/2 + t;
        elseif iDerivative == 2
            f = ones(size(t));
        elseif iDerivative > 2
            f = zeros(size(t));
        end
    else
        disp(['iShape = ', int2str(iShape)])
        disp('This program is not for this value of iShape')
        disp('Terminate this program by Ctrl+C')
        pause
    end
elseif iDegree == 3
    if iShape == 1
        if iDerivative == 0
            f = (1/16)*(-1+t+9*t.^2-9*t.^3);
        elseif iDerivative == 1
            f = (1/16)*(1+18*t-27*t.^2);
        elseif iDerivative == 2
            f = (9/8)*(1-3.*t);
        elseif iDerivative == 3
            f = -(27/8)*ones(size(t));
        elseif iDerivative > 3
            f = zeros(size(t));
        end
    elseif iShape == 2
        if iDerivative == 0
            f = (9/16)*(1-3*t-t.^2+3*t.^3);
        elseif iDerivative == 1
            f = (9/16)*(-3-2*t+9*t.^2);
        elseif iDerivative == 2
            f = (9/8)*(-1+9.*t);
        elseif iDerivative == 3
            f = (81/8)*ones(size(t));
        elseif iDerivative > 3
            f = zeros(size(t));
        end
    elseif iShape == 3
        if iDerivative == 0
            f = -(9/16)*(-1 - 3*t + t.^2 + 3*t.^3);
        elseif iDerivative == 1
            f = -(9/16)*(-3 + 2*t + 9*t.^2);
        elseif iDerivative == 2
            f = -(9/8)*(1+9*t);
        elseif iDerivative == 3
            f = -(81/8)*ones(size(t));
        elseif iDerivative > 3
            f = zeros(size(t));
        end
    elseif iShape == 4
        if iDerivative == 0
            f = (1/16)*(-1-t+9*t.^2+9*t.^3);
        elseif iDerivative == 1
            f = (1/16)*(-1+18*t+27*t.^2);
        elseif iDerivative == 2
            f = (9/8)*(1+3*t);
        elseif iDerivative == 3
            f = (27/8)*ones(size(t));
        elseif iDerivative > 3
            f = zeros(size(t));
        end
    else
        disp(['iShape = ', int2str(iShape)])
        disp('This program is not for this value of iShape')
        disp('Terminate this program by Ctrl+C')
        pause
    end
elseif iDegree == 4
    if iShape == 1
        if iDerivative == 0
            f = (1/6)*t.*(1 - t - 4*t.^2 + 4*t.^3);
        elseif iDerivative == 1
            f = (1/6)*(1 - 2*t - 12*t.^2 + 16*t.^3);
        elseif iDerivative == 2
            f = -1/3 -4*t + 8*t.^2;
        elseif iDerivative == 3
            f = -4 + 16*t;
        elseif iDerivative == 4
            f = 16*ones(size(t));
        elseif iDerivative > 4
            f = zeros(size(t));
        end
    elseif iShape == 2
        if iDerivative == 0
            f = -(4/3)*t.*(1 - 2*t - t.^2 + 2*t.^3);
        elseif iDerivative == 1
            f = -(4/3)*(1 - 4*t - 3*t.^2 + 8*t.^3);
        elseif iDerivative == 2
            f = (16/3) + 8*t - 32*t.^2;
        elseif iDerivative == 3
            f = 8 - 64*t;
        elseif iDerivative == 4
            f = -64*ones(size(t));
        elseif iDerivative > 4
            f = zeros(size(t));
        end
    elseif iShape == 3
        if iDerivative == 0
            f = 1 -  5*t.^2 + 4*t.^4;
        elseif iDerivative == 1
            f = 2*t.*(-5 + 8*t.^2);
        elseif iDerivative == 2
            f = -10 + 48*t.^2;
        elseif iDerivative == 3
            f = 96*t;
        elseif iDerivative == 4
            f = 96*ones(size(t));
        elseif iDerivative > 4
            f = zeros(size(t));
        end
    elseif iShape == 4
        if iDerivative == 0
            f = -(4/3)*t.*(-1 -2*t + t.^2 + 2*t.^3);
        elseif iDerivative == 1
            f = -(4/3)*(-1 -4*t + 3*t.^2 + 8*t.^3);
        elseif iDerivative == 2
            f = 16/3 - 8*t -32*t.^2;
        elseif iDerivative == 3
            f = -8*(1 + 8*t);
        elseif iDerivative == 4
            f = -64*ones(size(t));
        elseif iDerivative > 4
            f = zeros(size(t));
        end
    elseif iShape == 5
        if iDerivative == 0
            f = (1/6)*t.*(-1 -t +4*t.^2 + 4*t.^3);
        elseif iDerivative == 1
            f = (1/6)*(-1 -2*t + 12*t.^2 + 16*t.^3);
        elseif iDerivative == 2
            f = -(1/3) + 4*t + 8*t.^2;
        elseif iDerivative == 3
            f = 4 + 16*t;
        elseif iDerivative == 4
            f = 16;
        elseif iDerivative > 4
            f = zeros(size(t));
        end
    end 
elseif iDegree > 4
    disp(['degree = ', int2str(iDegree)])
    disp('This program is not for this value of degree')
    disp('Terminate this program by Ctrl+C')
    pause
end
end
% 
