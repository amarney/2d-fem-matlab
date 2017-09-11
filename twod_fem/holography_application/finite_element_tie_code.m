%% FINITE ELEMENT for 2D TIE-CODE
% solution: 
%           phi(x,y) = 
% coeffs:  
%           I_AA, I_BB, I_CC intensity images 
% pde:         
%           \nabla \cdot (I \nabla \phi) = -k I_z
% weak: 
%           - \int I \phi_x v_x - \int I \phi_y v_y = -\int k I_z v
% matrix:
%           (I_{xx} + I_{yy})u = k_0 + b_{I_{xx}} + b_{I_{yy}}
%
% Copyright 2017 Angelo Marney, Department of Mathematics, Virginia Tech
% For questions or concerns, contact me at: amarney@vt.edu

%% Set up image data: (courtesy of Hongbo Zhang)
negativefile = 'data/-1.png';
zerofile = 'data/0.png';
positivefile = 'data/1.png';
I161 = (im2double(imread(negativefile)));
I171 = (im2double(imread(zerofile)));
I181 = (im2double(imread(positivefile)));
A=  [90:180]; B = [1:91]; 
I16=I161(A,B); I17=I171(A,B); I18=I181(A,B);
IAA=I16;  % I - Z 
IBB=I17;  % I0
ICC=I18;  % I + Z

%% Set up parameter and spacing: 
lambda=633e-9;   
k0=(2*pi)/lambda;
dz=0.001; dx = 4.5e-6; dy = 4.5e-6; nPixels = 91; 
% note: dx, dy, and nPixels are redefined in matrix and vector assemblers.
%           need to fix this (pass dx, dy, nPixels straight into matrix and
%           vector assemblers.)

%% Matrices to be used as coefficients:
iName = IBB;
izName = (ICC - IAA)/(2*dz);  

%% Make a data structure for the mesh:
xMin = 0; xMax = nPixels*dx; yMin = 0; yMax = nPixels*dy; nx = 20; ny = 20;
Mesh = mesh_generator_2d_triangle(xMin, xMax, yMin, yMax, nx, ny);

%% Add edge data to mesh:
[edge, M_e, M_ne] = mesh_edge_generator(Mesh,1,1);
Mesh.edge = edge; Mesh.M_e = M_e; Mesh.M_ne = M_ne;

%% Make a data structure for finite element shape function access:
iDegree = 2; Fem = fem_generator_2d_lagrange(Mesh, iDegree);
Fem1 = Fem; Fem2 = Fem;

%% Prepare dof_u and ess bc data structures to handle boundary condition:
bName = @(x,y) x.*(x - xMax).*y.*(y - yMax);
bdName = @(x,y) x.*(x - xMax).*y.*(y-  yMax);
[dof_u, nt] = bc_array_generator_2d(Fem, bName, bdName);
g_D = @(x,y) zeros(size(x)); % BOUNDARY CONDITION !
u_e = zeros(length(Fem.point), 1);
I = find(nt == 0); u_e(I) = g_D(Fem.point(1, I), Fem.point(2, I));

%% Make the matrices:
nQuadraturePoint = 7;

xDerivative1 = 1; yDerivative1 = 0; xDerivative2 = 1; yDerivative2 = 0;
Ih_xx = matrix_2d_image_global(iName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
I_xx = Ih_xx(dof_u, dof_u);

xDerivative1 = 0; yDerivative1 = 1; xDerivative2 = 0; yDerivative2 = 1;
Ih_yy = matrix_2d_image_global(iName, Mesh, Fem1, ... 
    xDerivative1, yDerivative1, Fem2, xDerivative2, yDerivative2, nQuadraturePoint);
I_yy = Ih_yy(dof_u, dof_u);

%% prepare the -k_0 vector:
izName = (ICC - IAA)/(2*dz);  
xDerivative = 0; yDerivative = 0;
kh_0 = vector_2d_image_global(k0*izName, Mesh, Fem, xDerivative, ... 
    yDerivative, nQuadraturePoint);
k_0 = kh_0(dof_u);

%% prepare the b/c vectors associated with the matrices:
tmp = -Ih_xx*u_e; b_I_xx= tmp(dof_u);
tmp = -Ih_yy*u_e; b_I_yy= tmp(dof_u);

%% solve the system:
%  (I_xx + I_yy)*u_fe = k_0 + b_I_xx +  b_Iy_y0
u = (I_xx + I_yy)\(k_0 + b_I_xx +  b_I_yy);
u_fe = u_e; u_fe(dof_u) = u;

%% Do a quick plot of the solution:
figure(1); clf;
trisurf(Fem.T', Fem.point(1, :), Fem.point(2, :), u_fe, ...
'facecolor','interp', 'EdgeColor','none')
xlabel('x'), ylabel('y'), zlabel('\phi (x,y,z_0)')


