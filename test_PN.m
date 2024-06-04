% Run test for even parity PN method for solving the radiative transfer equation 
% u_t + s.nabla_r u + sigma_t*u = sigma_s*int_Stheta(s.s')u(s')ds' + q
% on a cylindrical domain in 3D

% Physical parameters
r = 0.1;     % cylinder radius 
H = 0.5;     % cylinder length
mus = 1;     % scattering coefficient
mua = 0.1;   % absorption coefficient

% Computational parameters
N=3;         % maximum order of harmonics in numerical solution
tol = 1e-3;  % tolerance for the iterative solver
g = 0;       % anisotropy (assumes Henyey-Greenstein kernel)
h = 0.01;    % maximal mesh size

% Create geometry
gm = multicylinder(r,H);                                % create cylindrical model
model = createpde;                                      % create a generic PDE model
model.Geometry = gm;                                    % assign the geometry to a generic PDE model
generateMesh(model,'GeometricOrder','linear','Hmax',h); % create a mesh of linear elements
p=model.Mesh.Nodes;                                     % list of coordinates of the nodes
e=surftri(model.Mesh.Nodes',model.Mesh.Elements')';     % list of boundary elements
t=model.Mesh.Elements;                                  % list of elements

%----- plot geometry and mesh if desired -----
% figure
% pdegplot(model,'FaceLabels','on')
%Generate mesh and plot the mesh.

% figure
% pdemesh(model);
% --------------------------------------------

%----- possible refinement procedure for the mesh -----
% for i=1:maxlevel
%     [p,e,t]=refine3(p,e,t);
% end

%-------------------------------------------------------

fprintf('Geometry was setup in %d sec\n',toc);

% obtain sizes of the problem
nt=size(t,2);                           % number of spatial triangles
np=size(p,2);                           % number of spatial points
Ne = sum(2*(0:2:N)+1);                  % number of even harmonics
No = sum(2*(1:2:N)+1);                  % number of odd harmonics
x = p(1,:); y = p(2,:); z = p(3,:);     % coordinates of the nodes

fprintf('Ne=%d, No=%d, N=%d, np=%d, nt=%d, dofs %d\n',Ne,No,N,np,nt,sum(2*(0:2:N)+1)*np)

% -------------- Angular assembly ------------------
 
fprintf('Assembling angular part ...\n');
tic

[A1,A2,A3] = getAngularA(N);   % assemble the angular part of the transport erm

% assemble the scattering integral using a legendre expansion of the phase
% function and the addition theorem for spherical harmonics 
thetaEven = zeros(Ne,1);
thetaOdd = zeros(No,1);
iter=1;
for l=0:2:N
    for m=-l:l
        thetaEven(iter) = g^l;
        iter=iter+1;
    end
end

iter=1;
for l=1:2:N
    for m=-l:l
        thetaOdd(iter) = g^l;
        iter=iter+1;
    end
end

%assemble full matrices for the scattering integral and absorption part
T1 = diag(sparse(thetaEven)); T2 = diag(sparse(thetaOdd)); 
I1 = speye(Ne); I2 = speye(No);


fprintf('Finished angular assembly in %d sec\n',toc)

%------------- Spatial assembly -------------------

fprintf('Assembling spatial part ...\n');
tic

[dx,dy,D3,detB]=assemble_A_dxdydz(p,t); %assemble spatial part of transport term

% Combine x and y components for convenience
D1=(dx+1i*dy);  
D2=(dx-1i*dy);

%assemble spacial mass matrices and stiffness matrices
[Laplace,M1]=assem3D(p,t,1,1);
[~,Mua1]=assem3D(p,t,[],mua);
[~,Mus1]=assem3D(p,t,[],mus);

volTet = detB/6; %volume of all the tetrahedra

% calculate inverse of odd scattering matrices
C2inv = getC2Inverse(N,mua,mus,thetaOdd,volTet);

fprintf('Finished spatial assembly in %d sec\n',toc)

% --------- Boundary assembley ---------
fprintf('Assembling boundary matrix ...\n');
tic
B_M = assemble_boundary(p,e,100,N);
R=assemble_robin2d(p,e);
fprintf('Finished boundary assembly in %d sec\n',toc)


%assemble source term
Q = zeros(np,Ne);
Q(:,1) = M1*exp(-( (x).^2+(y).^2) * 20)';

%Define system to be solved
A1s = A1';A2s = A2';A3s = A3';
D1s = D1';D2s = D2';D3s = D3';

A = @(x) D1 * x * A1s + D2 * x * A2s + D3 * x * A3s;
As= @(x) D1s* x * A1  + D2s * x * A2  + D3s * x * A3;

C1 = @(x) Mua1 * x + Mus1 * x * (I1-T1)';
C2i=@(x) C2inv.*x;
K = @(x) As(C2i(A(x))) + C1(x) + reshape(B_M*x(:),np,Ne);

% -------- Preconditioning -------------------
fprintf('Factoring preconditioner ...\n');
tic
[L,U,P,Qlu] = lu(Laplace+M1+R);
fprintf('Factoring succesfull in %d sec\n',toc)


% --------- Solving the system ----------------
fprintf('Solving ....\n');tic;
pre=@(x) Qlu*(U\(L\(P*x)));

[phi,flag,relres,it]=mypcg(K,Q,tol,500,pre);
if flag==0
        fprintf('pcg converged in %f sec and %d it to desired tol %e with relres=%e\n',toc,it,tol,relres);
else
    fprintf('pcg did not converge to desired tol in %f sec and %d it with relres=%e\n',toc,it,relres);
end

