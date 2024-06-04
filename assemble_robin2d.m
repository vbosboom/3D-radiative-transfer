function [robin,proj_lin2const,mass_pcw_const]=assemble_robin2d(p,e)
%proj_lin2const int phi_i chi_k dGamma projects piecwise linear functions in the domain to piecewise constant on the boundary.
%robin          int phi_i phi_k dGAmma mass matrix of basis functions at the boundary
%mass_pcw_const int chi_i chi_k dGamma mass matrix of piecwise constants on the boundary

% The Basisfunctions for the reference Intervall:
ne=size(e,2);   %number of edge surfaces
np=size(p,2);   %number of points 

Psi=@(x) [1-x(1)-x(2);x(1);x(2)];
% dPsi=[-1,1,0;-1,0,1];
%Integrate all the basis functions and their derivatives over the Unittriangle 

no_basis = 3;

precision = 2;
[ ip, w ] = integrationrule2D( precision );

%Integrate all the basis functions and their derivatives over the reference tetrahedron
Psi_eval = zeros( no_basis, length( w ) );
for i = 1:length( w )
    Psi_eval( :, i ) = Psi( ip( i, : ) );
end

int_psi_prod = Psi_eval * ( repmat( w, 1, no_basis ) .* Psi_eval' );
int_psi = Psi_eval * w;

p1=p(1:3,e(1,:)); %first node of faces
p2=p(1:3,e(2,:)); %second node of faces
p3=p(1:3,e(3,:)); %third node of faces

normal=cross(p3-p1,p2-p1);
detB=sqrt(sum(normal.^2,1));

row=kron(e(1:3,:),ones(1,3));
col=reshape(repmat(reshape(e(1:3,:),1,3*ne),3,1),3,3*ne);
val=kron(detB,int_psi_prod);
robin=sparse(row,col,val,np,np);

mass_pcw_const=diag(sparse(abs(detB)))/2;    %int_T chi_i chi_k dGamma =1/2 for T reference triangle.

row=repmat(1:ne,3,1);
col=e(1:3,:);
val=kron(detB,int_psi);
proj_lin2const=sparse(row,col,val,ne,np);

end

