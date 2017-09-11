function f = shape_function_generator_2d_lagrange_triangle_reference(hx, hy, ...
    iDegree, hxDerivative, hyDerivative, iShape)
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu
 if iDegree == 0
     if iShape == 1
         if iDerivative == 0
             f = ones(size(hx));
         elseif iDerivative > 0
             f = zeros(size(hx));
         end
     else
         disp(['iShape = ', int2str(iShape)])
         disp('This program is not for this value of iShape')
         disp('Terminate this program by Ctrl+C')
         pause
     end        
 elseif iDegree == 1
    if iShape == 1
        if hxDerivative == 0 && hyDerivative == 0
            f = 1 - hx - hy;
        elseif hxDerivative == 1 && hyDerivative == 0
            f = -1*ones(size(hx));
        elseif hxDerivative == 0 && hyDerivative == 1
            f = -1*ones(size(hx));
        elseif hxDerivative > 1 || hyDerivative > 1
            f = zeros(size(hx));
        end
    elseif iShape == 2
        if hxDerivative == 0 && hyDerivative == 0
            f = hx;
        elseif hxDerivative == 1 && hyDerivative == 0
            f = 1*ones(size(hx));
        elseif hxDerivative == 0 && hyDerivative == 1
            f = zeros(size(hx));
        elseif hxDerivative > 1 || hyDerivative > 1
            f = zeros(size(hx));
        end
    elseif iShape == 3
        if hxDerivative == 0 && hyDerivative == 0
            f = hy;
        elseif hxDerivative == 1 && hyDerivative == 0
            f = zeros(size(hx));
        elseif hxDerivative == 0 && hyDerivative == 1
            f = 1*ones(size(hx));
        elseif hxDerivative > 1 || hyDerivative > 1
            f = zeros(size(hx));
        end
    else
        disp(['iShape = ', int2str(iShape)])
        disp('This program is not for this value of iShape')
        disp('Terminate this program by Ctrl+C')
        pause
    end           
 elseif iDegree == 2
     if iShape == 1
         if hxDerivative == 0 && hyDerivative == 0
             f = 2*(1 - hx - hy).*((1 - hx - hy) - .5);
         elseif hxDerivative == 1 && hyDerivative == 0
             f = -3*ones(size(hx)) + 4*hy + 4*hx;
         elseif hxDerivative == 0 && hyDerivative == 1
             f = -3*ones(size(hx)) + 4*hx + 4*hy;
         elseif hxDerivative == 1 && hyDerivative == 1
             f = 4*ones(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));
         end
     elseif iShape == 2
         if hxDerivative == 0 && hyDerivative == 0
             f = 2*hx.*(hx - .5);
         elseif hxDerivative == 1 && hyDerivative == 0
             f = 4*hx - ones(size(hx));
         elseif hxDerivative == 0 && hyDerivative == 1
             f = zeros(size(hx));
         elseif hxDerivative == 1 && hyDerivative == 1
             f = zeros(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));
         end
     elseif iShape == 3
         if hxDerivative == 0 && hyDerivative == 0
             f = 2*hy.*(hy - .5);
         elseif hxDerivative == 1 && hyDerivative == 0
             f = zeros(size(hx));
         elseif hxDerivative == 0 && hyDerivative == 1
             f = 4*hy - ones(size(hx));
         elseif hxDerivative == 1 && hyDerivative == 1
             f = zeros(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));             
         end
     elseif iShape == 4
         if hxDerivative == 0 && hyDerivative == 0
             f = 4*(1 - hx - hy).*hx;
         elseif hxDerivative == 1 && hyDerivative == 0
             f = 4*ones(size(hx)) - 8*hx - 4*hy;
         elseif hxDerivative == 0 && hyDerivative == 1
             f = -4*hx;
         elseif hxDerivative == 1 && hyDerivative == 1
             f = -4*ones(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));              
         end
     elseif iShape == 5
         if hxDerivative  ==0 && hyDerivative == 0
             f = 4*hy.*hx;
         elseif hxDerivative == 1 && hyDerivative == 0
             f = 4*hy;
         elseif hxDerivative == 0 && hyDerivative == 1
             f = 4*hx;
         elseif hxDerivative == 1 && hyDerivative == 1
             f = 4*ones(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));   
         end
     elseif iShape == 6
         if hxDerivative == 0 && hyDerivative == 0
             f = 4*hy.*(1 - hx - hy);
         elseif hxDerivative == 1 && hyDerivative == 0
             f = -4*hy;
         elseif hxDerivative == 0 && hyDerivative == 1
             f = 4*ones(size(hx)) - 8*hy - 4*hx;
         elseif hxDerivative == 1 && hyDerivative == 1
             f = -4*ones(size(hx));
         elseif hxDerivative > 1 || hyDerivative > 1
             f = zeros(size(hx));                
         end
        else
            disp(['iShape = ', int2str(iShape)])
            disp('This program is not for this value of iShape')
            disp('Terminate this program by Ctrl+C')
            pause
     end   
end