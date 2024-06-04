%compute the stiffness and mass matrix over the spatial grid
%ks stiffness matrix int_R grad diffusion * phi_i * grad phi_j dR
%km mass matrix int_R absorption * phi_i phi_j dR
function [ks,km]=assem3D(p,t,diffusion,absorption)


nt=size(t,2);   %number of tetrahedrons
np=size(p,2);   %number of nodes

% if no coefficients are mentioned
if isempty(diffusion) && isempty(absorption)
    ks = sparse(np,np);
    km = sparse(np,np);
    return
end

% The basis functions on the reference element
Psi=@(x) [1-x(1)-x(2)-x(3); x(1); x(2);x(3)]; 
dPsi=[-1,1,0,0;
      -1,0,1,0;
      -1,0,0,1];
no_basis=4;
precision=2;
[ip,w]=intrule_tet(precision);

%Integrate all the basis functions and their derivatives over the reference tetrahedron 
Psi_eval=zeros(no_basis,length(w));
for i=1:length(w)
    Psi_eval(:,i)=Psi(ip(i,:));
end
    
int_Psi_prod=Psi_eval*(repmat(w,1,no_basis).*Psi_eval');

row_pos=kron(t(1:4,:),ones(1,4));
col_pos=repmat(reshape(t(1:4,:),1,4*nt),4,1);

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
detB=sum(a.*bxc);

if isempty(absorption)
    km = sparse(np,np);
else
    km=sparse(row_pos,col_pos, kron(absorption.*abs(detB),int_Psi_prod),np,np);
end

%B_inverse=[bxc'.*v;axc'.*v;axb'.*v]; v=[1 -1 1];
B_inverse=reshape([ bxc(1,:)./detB; -axc(1,:)./detB; axb(1,:)./detB; ...
                    bxc(2,:)./detB; -axc(2,:)./detB; axb(2,:)./detB; ...
                    bxc(3,:)./detB; -axc(3,:)./detB; axb(3,:)./detB] ...
                    , 3,3*nt);


if isempty(diffusion)
    ks = sparse(np,np);
else
    BIdPsi=B_inverse'*dPsi;

    BIdPsi = reshape( BIdPsi, [ 3 nt 4] );
    BIdPsi = permute( BIdPsi, [ 1 3 2 ] );
    BIdPsi = reshape( BIdPsi, [ 3 4*nt ] );

    %calculate (B^-T grad)^T psi B^-T grad psi
    a=BIdPsi(:,1:4:end);
    b=BIdPsi(:,2:4:end);
    c=BIdPsi(:,3:4:end);
    d=BIdPsi(:,4:4:end);
    AdetB=(abs(detB).*diffusion)/6; %1/6 mass of unit triangle
    sp1=sum(a.*a,1).*AdetB;
    sp2=sum(a.*b,1).*AdetB;
    sp3=sum(a.*c,1).*AdetB;
    sp4=sum(a.*d,1).*AdetB;
    sp5=sum(b.*b,1).*AdetB;
    sp6=sum(b.*c,1).*AdetB;
    sp7=sum(b.*d,1).*AdetB;
    sp8=sum(c.*c,1).*AdetB;
    sp9=sum(c.*d,1).*AdetB;
    sp10=sum(d.*d,1).*AdetB;
    ks=reshape([sp1;sp2;sp3;sp4;sp2;sp5;sp6;sp7;sp3;sp6;sp8;sp9;sp4;sp7;sp9;sp10], 4, 4*nt);
    ks=sparse(row_pos,col_pos,ks,np,np);
end