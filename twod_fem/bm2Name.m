function bm2Name = bm2Name(x,y,xInterpolate,wInterpolate)
beps = 10e-7;
if ((x >= 0) & (x <= 1) & (y == 1)) % top
    bm2Name = 0;
elseif ((x == 0) & (y >= 0) & (y <= 1)) % left
    bm2Name = 0;
elseif ((x == 1) & (y >= 0) & (y <= 1)) % right
    bm2Name = 0;
else % else, on bottom
    %bm2Name = interp1(xInterpolate, wInterpolate, x, 'pchip');
    bm2Name = interp2(xInterpolate,wInterpolate,wInterpolate, x,y,'cubic')
end
% g2_D = @(x,y) 0.*((x >= 0) & (x <= 1) & (y == 1)) + ...
%     0.*((x == 0) & (y >= 0) & (y <= 1)) + ...
%     0.*((x == 1) & (y >= 0) & (y <= 1)) + ...
%     interp1(xInterpolate, wInterpolate, x, 'pchip').*((x >= 0) & (x <= 1) & (y <= );
end
