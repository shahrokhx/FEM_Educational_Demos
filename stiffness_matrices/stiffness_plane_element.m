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
% Calculating a single 4-nodes plane stress element stiffness matrix by 
% numerical integration
%
% (Logan - Example 10.4 - page 511)
%
%        |           
%        |
%        |  (3,4) O------O (5,4)
%        |        |      |
%        |        |      |
%        |  (3,2) O------O (5,2)
%        |
%        +--------------------------
%
%
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
nElemNode = 4; % number of nodes per elements

% point1 = (3,2)  point2 = (5,2)  point3 = (5,4)  point4 = (3,4)
%              [1]   [2]   [3]   [4]
coordinates = [ 3     5     5     3     %  x
                2     2     4     4]    %  y
            

% material properties             
E  = 30e6;
nu = 0.25;
h  = 1;

%% Elasticity Matrix
D = (E/(1-nu^2))*    [1   nu    0     ; 
                      nu  1     0     ; 
                      0   0  (1-nu)/2]
                  
%% Numerical Integration Setup
xiList = [-0.57735,    0.57735,   -0.57735,    0.57735  ;
          -0.57735,   -0.57735,    0.57735,    0.57735] ;
      
wList  = [1.0     ,    1.0    ,    1.0    ,    1.0    ] ;
                  
%% Stiffness Matrix Calculation by Numerical Integration
K = zeros(nElemNode * nDof);

for intpt = 1 : length(xiList)
  fprintf('----------------------------------------------------\n')
  xi    = xiList(:,intpt)';  % local corrdinate
  w     = wList (intpt)   ;
  
  fprintf('[%d]    xi = %f     w = %g \n',intpt,xi,w)
  % shape functions and derivatives at local integration points
  N     = shapefunctions     (nDim,nElemNode,xi) 
  dNdxi = shapefunctionderivs(nDim,nElemNode,xi)
  
  % Jacobian matrix
  J     = coordinates * dNdxi
  
  % inverse and determinent of the Jacobian matrix
  invJ  = inv(J)
  detJ  = abs(det(J))
  
  % derivatives of the shape functions in the global coordinate 
  dNdx  = dNdxi * invJ
    
  % B matrix
  B = zeros(3, nElemNode*nDof);
  for i = 1 : nElemNode
      idx = 2*(i-1) + 1;
      B(1,idx  ) = dNdx(i,1);
      B(3,idx+1) = dNdx(i,1);
      B(2,idx+1) = dNdx(i,2);
      B(3,idx  ) = dNdx(i,2);
  end
  
  % numerical summation
  K     = K     +        B'*D*B  *  detJ  *  w;
  
end
fprintf('----------------------------------------------------\n')
K