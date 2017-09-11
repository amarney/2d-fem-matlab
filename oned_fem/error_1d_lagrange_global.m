function errorGlobal = error_1d_lagrange_global(uName, uGlobal, Mesh, Fem, iDerivative, nGaussPoint)
iDegree = Fem.degree;
%--------------------------------------------------------------------------
%error_1d_lagrange_global:
%   computes the global error on a given domain using the L^2 norm,
%
%   norm(u - u_h)_{L^2(\omega) =  \sum (\int_{K} (u^{(r)}(x) - u_h^{(r)}(x))^2) dx.   
%
%   Can be used to numerically identify the order of convergence for the
%   finite element interpolation. See example below.
%
%   Parameters  :   uName - function handle for true function.
%                   uGlobal - array containing coefficients of finite
%                             element solution (matches solution at nodes).
%                   Mesh - the Mesh data structure which contains node
%                          coordinates, the element connectivity array,
%                          and the edge connectivity array.
%                   Fem - the Fem data structure which contains
%                         x-coordinates of of the finite element space, the
%                         connectivity array, and the degree of the basis
%                         polynomials.           
%                   iDerivative - index specifying which derivative of the
%                                 shape function is evaluated.                   
%                   nGaussPoint - integer representing the number of gauss
%                                  to be generator.
%
%   Return      :   errorGLobal - integer representing the global error.
%   
%   Example     :
%       To find the global error of u(x) = cos(x) over a domain [0,1] using
%       second degree basis polynomials with 5 elements do
%
%       domain = [0,1]; nElement = 5;
%       Mesh = mesh_generator_1d(domain, nElement);
%       iDegree = 2; Fem = fem_generator_1d_lagrange(Mesh, iDegree);
%       uName = @(x) cos(x); uGlobal = cos(Fem.point); 
%       iDerivative = 0; nGaussPoint = 5;      
%       errorGlobal = error_1d_lagrange_global(uName, uGlobal, Mesh, Fem, ... 
%           iDerivative, nGaussPoint) 
%               
%       To numerically identify the order of convergence for the finite
%       element interpolation using second degree basis polynomials over
%       the domain [0,1], do
%   
%       uName = @(x) cos(x); uPrimeName = @(x) -sin(x);
%       domain = [0,1]; iDegree = 2; nGaussPoint = 5;
%       nElements = [10,20,30,40,50,60];
%       errorL2 = zeros(size(nElements)); 
%       errorH1 = zeros(size(nElements));
%       for i = 1:length(nElements)
%           nElement = nElements(i); 
%           Mesh = mesh_generator_1d(domain, nElement);
%           Fem = fem_generator_1d_lagrange(Mesh, iDegree);
%           uGlobal = feval(uName,Fem.point);
%           errorL2(i) = sqrt(error_1d_lagrange_global(uName, uGlobal, Mesh, ...
%               Fem, 0, nGaussPoint));
%           errorH1(i) = sqrt(error_1d_lagrange_global(uPrimeName, uGlobal, Mesh, ...
%               Fem, 1, nGaussPoint));
%       end
%       h = 1./nElements %spacing
%       rc = polyfit(log(h), log(errorL2),1);
%       r_L2 = rc(1), C_L2 = rc(2)              C_L2 = exp(rc(2))?
%       rc = polyfit(log(h), log(errorH1),1);
%       r_H1 = rc(1), C_H1 = rc(2)              C_H1 = exp(rc(2))?
%
%       % here r is the order of convergence
%       % errorL2(h) = C*h^r
%       % errorH1(h) = C*h^r
%--------------------------------------------------------------------------
errorGlobal = 0;
for k = 1:size(Mesh.element,2)
    element = Mesh.node(Mesh.element(:,k));
    uLocal = uGlobal(Fem.T(:,k));
    errorLocal = error_1d_lagrange_local(uName, uLocal, ...
        element, iDegree, iDerivative, nGaussPoint);
    errorGlobal = errorGlobal + errorLocal;
end
return;