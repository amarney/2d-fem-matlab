function point_is_in = ptInTriangle(pt, element)
    p1.x = element(1,1); p1.y = element(2,1);
    p2.x = element(1,2); p2.y = element(2,2);
    p3.x = element(1,3); p3.y = element(2,3);
    A = 1/2 * (-p2.y .* p3.x + p1.y .* (-p2.x + p3.x) + p1.x .* (p2.y - p3.y) + p2.x .* p3.y);
    if (A < 0)
        signt = -1;
    else
        signt = 1;
    end
    s = (p1.y .* p3.x - p1.x .* p3.y + (p3.y - p1.y) .* pt.x + (p1.x - p3.x) .* pt.y) .* signt;
    t = (p1.x .* p2.y - p1.y .* p2.x + (p1.y - p2.y) .* pt.x + (p2.x - p1.x) .* pt.y) .* signt;
    point_is_in = s >= 0 & t >= 0 & (s + t) <= 2 .* A .* signt;
end
% added .* instead of * to handle array of points
% aded & instead of && to handle array of points
%http://jsfiddle.net/PerroAZUL/zdaY8/1/