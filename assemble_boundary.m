function B_Me = assemble_boundary(p,e,precision,N)

Ne = sum(2*(0:2:N)+1); %number of even moments
np = size(p,2); %number of points
nBt = size(e,2); %number of boundary elements

%Gauss-Legendre integration points and weights in 1D
[points, weights] = lgwt(precision,-1,1);

%Initialize empty matrices
B_Me = sparse(np*Ne,np*Ne); %full matrix
Phi_Mat = zeros(nBt*3,3); %spatial mass matrix

Deltabnd = getBoundaryAreas(p,e); %result is twice the area of the boundary triangle

%Construct spatial mass matrix for all elements
for ind1 = 1:3
    for ind2 = 1:3
        MBelem = Deltabnd./24.*(1+(ind1==ind2));
        Phi_Mat(ind1:3:end,ind2) = MBelem;
    end
end

%Indices for spatial mass matrix
row_Phi = reshape(kron(ones(1,3),e),[],3);
col_Phi = kron(e',ones(3,1));

x = p(1,:)';
y = p(2,:)';
z = p(3,:)';

%calculate the outward normal and the angular coordinates
p1 = [x(e(1,:)),y(e(1,:)),z(e(1,:))];
p2 = [x(e(2,:)),y(e(2,:)),z(e(2,:))];
p3 = [x(e(3,:)),y(e(3,:)),z(e(3,:))];

normal = cross(p2-p1,p3-p1);
detB=sqrt(sum(normal.^2,2)); 
n=normal./repmat(detB,1,3); 

n = n.*(1-2.*(sum((p1+n)>2 | (p1+n)<0,2)==0)); %CHECK IF NORMAL IS ORIENTED THE RIGHT WAY

cos_t_n = n(:,3);
sin_t_n = sqrt(n(:,1).^2+n(:,2).^2);
phi_n=atan2(n(:,1),n(:,2)); %angle of the norm

%functions for in integrals
G1 = @(B) 2.*cos(B).*sin(B)+2.*B-pi; %for the case k^2=1
G2 = @(k,B) 2./(k-1).*sin((k-1).*B)+2./(k+1).*sin((k+1).*B); %for the case k^2 ~=1

F1 = @(B) 4.*B-2.*pi; %for the case k=0
F2 = @(k,B) 4./k.*sin(B.*k); %for the case k~=0

alpha = sqrt(1-points.^2).'.*sin_t_n;
beta = points.'.*cos_t_n;
gamma = -beta./alpha;

B = zeros(size(alpha));

B(abs(beta)<alpha) = acos(gamma(abs(beta)<alpha));
B(beta>alpha) = pi;

%obtain angular integral
%loop over even harmonics
for l= 0:2:N
    P1 = legendre(l,points);
    for k=0:2:N
        P2 = legendre(k,points);
        for m=-l:l
            Plm=P1(abs(m)+1,:); %get correct P_l^m
            for n=-k:k
                Pkn=P2(abs(n)+1,:); %get correct P_k^n

                f1=sqrt((2*l+1)/(4*pi) * factorial(l-abs(m))/factorial(l+abs(m))) * (-1)^(0.5*(m+abs(m))); %prefactors
                f2=sqrt((2*k+1)/(4*pi) * factorial(k-abs(n))/factorial(k+abs(n))) * (-1)^(0.5*(n+abs(n)));
                
                %indices for the spherical harmonics in combined matrix
                col_off = (momentToDof(l,m)-1)*np;
                row_off = (momentToDof(k,n)-1)*np;

                int1=0;
                if (m==n)
                    int1 = int1+beta.*F1(B);
                else
                    int1 = int1+beta.*F2(m-n,B);
                end

                if (abs(m-n)==1)
                    int1 = int1 + alpha.*G1(B);
                else
                    int1 = int1 + alpha.*G2(m-n,B); %angular integral
                end

                val=int1.*f1.*f2.*exp(1i.*(n-m).*phi_n).*Plm.*Pkn*weights;
                val = kron(val,ones(3)).*Phi_Mat;

                Col = col_off+col_Phi;
                Row = row_off+row_Phi;

                
                B_Me = B_Me+sparse(Row,Col,val,np*Ne,np*Ne); %collect terms in full matrix
            end
        end
    end
end