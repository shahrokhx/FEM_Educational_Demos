%%               E D U C A T I O N A L      S N I P P E T S                     
%__________________________________________________________________________
% 
%                       Finite Element Methods
%                     Developed by SHAHROKH SHAHI 
%                           (www.sshahi.com)
%
%                   Georgia Institute of Technology
%__________________________________________________________________________
%
%{ 
This simple script demonstrates how the stiffness matrix of a "truss" 
element can be obtained by numerical integration.

%}
%% Initial Setup 
clc
clear
close all

format short g
format compact

% adding library folder(s) to the path
path = mfilename('fullpath');
path(end-length(mfilename):end)=[];
addpath(fullfile(path,'lib_integration'))
%% General Information

nDim = 2;      % number of dimensions (2D)
nDof = 2;      % number of degrees of freedom
nElemNode = 2; % number of nodes per elements

% point1 = (4,0)  point2 = (0,3)
%               [1]    [2]
coordinates = [  4      0    ;   %  x
                 0      3   ]    %  y

% material properties             
E = 210e9
A = 6e-4


%% Numerical Integration Setup
xiList = [-0.5773502692,  0.5773502692];
wList  = [1.0          ,  1.0         ];

xiList = 0;
wList  = 2.0;

%% Numerical Integration Procedure

dCoord = diff(coordinates');
L      = norm(dCoord);
coordinates = [0 L]; % in 1D view
c = dCoord(1)/L;
s = dCoord(2)/L;


k = zeros(nDof);
for intpt = 1 : length(xiList)
  fprintf('----------------------------------------------------\n')
  xi    = xiList(intpt);  % local corrdinate
  w     = wList (intpt);
  
  fprintf('[%d]    xi = %f     w = %g \n',intpt,xi,w)
  % shape functions and derivatives at local integration points
  N     = shapefunctions     (nDim-1,nElemNode,xi) 
  dNdxi = shapefunctionderivs(nDim-1,nElemNode,xi)
  
  % Jacobian matrix
  J     = coordinates * dNdxi
  
  % inverse and determinent of the Jacobian matrix
  invJ  = inv(J)
  detJ  = abs(det(J))
  
  % derivatives of the shape functions in the global coordinate 
  dNdx  = dNdxi * invJ
    
  % B matrix
  B     = dNdx'
  
  % numerical summation
  k     = k     +        B'*E*B  *  detJ * A * w;
  
end
fprintf('----------------------------------------------------\n')
k


fprintf('----------------------------------------------------\n')
fprintf('In Global Coordinates: \n')

% 2-D coordinates
k4 = zeros(nElemNode * nDof);
k4([1,3],[1,3]) = k;
T = [[c s; -s c] zeros(2)    ; 
     zeros(2)    [c s; -s c]];
K = T' * k4 * T
fprintf('----------------------------------------------------\n')

