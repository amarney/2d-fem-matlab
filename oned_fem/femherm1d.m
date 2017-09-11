function Fem = femherm1d(Mesh, iDegree)
hermDegree = iDegree - 2; % iDegree = 3 is 1st order hermite space
nElement = size(Mesh.element,2);
t = zeros(hermDegree+1,nElement); %t(:,1) = (1:hermDegree+1)';
t(:,1) = [1, hermDegree + 1]';
for k=2:nElement
    t(:,k) = t(:,k-1)+hermDegree;
end
%t(:,(nElement/2 + 1):nElement) = t(:,1:nElement/2 );
p = zeros(1, (hermDegree+1) + (nElement - 1)*hermDegree);
k = 1;
element = Mesh.node(Mesh.element(:,k));
h = (element(2) - element(1))/hermDegree;
p(1:hermDegree+1) = element(1):h:element(2);
pointCount = hermDegree + 1;
for k=2:nElement
    element = Mesh.node(Mesh.element(:,k));
    h = (element(2) - element(1))/hermDegree;
    p(pointCount+1:pointCount+hermDegree) = element(1)+h:h:element(2);
    pointCount = pointCount+hermDegree;
end
p = [p,p];
t = [t; t+nElement+1];
%t(1,:) = t(1,:) + nElement;
%t(3,:) = t(3,:) + nElement;
%foo = find(t>42); t(foo) = t(foo) - nElement;
%t = [t;t + length(Mesh.node)];
Fem = struct('point', p, 'T', t, 'degree', hermDegree);
end
