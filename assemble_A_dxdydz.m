%Assembly routine of the transport matrices on the spatial domain
%dx = int_R chi_i * dphi_j/dx dR
%dy = int_R chi_i * dphi_j/dy dR
%dz = int_R chi_i * dphi_j/dz dR

% here chi_i is a piecewise constant basis function and phi_i a piecewise
% linear basis function
function [dx,dy,dz,detB]=assemble_A_dxdydz(p,t)


nt=size(t,2);   %number of tetrahedrons
np=size(p,2);   %number of points

% gradients of the basis functions on the reference interval
dPsi=[-1,1,0,0;
      -1,0,1,0;
      -1,0,0,1];

 
%calculate volumes of tetrahedrons
p1=p(1:3,t(1,:));
p2=p(1:3,t(2,:));
p3=p(1:3,t(3,:));
p4=p(1:3,t(4,:));
a=p2-p1;
b=p3-p1;
c=p4-p1;

bxc=cross(b,c);
axc=cross(a,c);
axb=cross(a,b);
detB=sum(a.*bxc);  % this is 6 times the volume of the tetrahedra

%mapping back to the reference element
B_inverse=reshape([ bxc(1,:)./detB; -axc(1,:)./detB; axb(1,:)./detB; ...
                    bxc(2,:)./detB; -axc(2,:)./detB; axb(2,:)./detB; ...
                    bxc(3,:)./detB; -axc(3,:)./detB; axb(3,:)./detB] ...
                    , 3,3*nt);

BIdPsi = B_inverse'*dPsi;
BIdPsi = reshape( BIdPsi, [ 3 nt 4] );
BIdPsi = permute( BIdPsi, [ 1 3 2 ] );
BIdPsi = reshape( BIdPsi, [ 3 4*nt ] );

% calculate the integral as int_t B^(-T) Grad Psi
detB_rep=reshape(repmat(abs(detB),4,1),1,4*nt);
dx = (BIdPsi(1,:).*detB_rep)/6;
dy = (BIdPsi(2,:).*detB_rep)/6;
dz = (BIdPsi(3,:).*detB_rep)/6;

row=repmat(1:nt,4,1);
col=t(1:4,:);

%assemble the final matrices
dx=sparse(row,col,dx,nt,np);
dy=sparse(row,col,dy,nt,np);
dz=sparse(row,col,dz,nt,np);

end