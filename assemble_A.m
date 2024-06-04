function A=assemble_A(N,mesh)

p=mesh.p;
t=mesh.t;
np=size(p,2);
nt=size(t,2);
[dx,dy,dz]=assemble_A_dxdydz(p,t);
unknowns=(N+1)^2*nnz(dx)*8;
I=sqrt(-1);
row_pos=zeros(unknowns,1);
col_pos=zeros(unknowns,1);
values=zeros(unknowns,1);
D1=(dx+I*dy);
D2=(dx-I*dy);
ind=0;
for l=1:2:N
    for m=-l+1:l-1
        off_row = l*(l-1)/2*nt + (m+l-1)*nt;
        off_col = (l-1)*(l-2)/2*np + (m+l-1)*np;
%         [l,1-m; l,m;  l,1+m]
        [i,j,val]=find([-a(l,1-m)*D1;  b(l,m)*dz; a(l,1+m)*D2]);
        nn=length(i);
        row_pos(ind+1:ind+nn)=off_row+i;
        col_pos(ind+1:ind+nn)=off_col+j;
        values(ind+1:ind+nn)=val;
        ind=ind+nn;
    end  
end
for l=1:2:N-1
     for m=-l:l
        off_row = l*(l-1)/2*nt + (m+l)*nt;
        off_col = l*(l+1)/2*np + (m+l)*np;
        if l==3, [l+1,1-m,l+1,m,l+1,m+1];, end
        [i,j,val]=find([-a(l+1,1-m)*D2,  b(l+1,m)*dz, a(l+1,m+1)*D1]);
        nn=length(i);        
        row_pos(ind+1:ind+nn)=off_row+i;
        col_pos(ind+1:ind+nn)=off_col+j;
        values(ind+1:ind+nn)=val;
        ind=ind+nn;
     end
end

k1=find(row_pos,1,'last');
k2=find(col_pos,1,'last');
k3=find(values,1,'last');
if(k1>unknowns), fprintf('Warning: allocate more memory for creating A\n');end;
if (k1==k2 && k2==k3)
    row_pos=row_pos(1:k1);
    col_pos=col_pos(1:k2);
    values=values(1:k3);
    A=sparse(row_pos,col_pos,values,nt*(N+2)*(N+1)/2, np*N*(N+1)/2, k3);
else
    fprintf('Error in Creating Differential operator\n');
end
end

function erg=a(l,m)
erg=0.5*((l+m)*(l+m-1)/((2*l+1)*(2*l-1)))^0.5;
end

function erg=b(l,m)
erg= ( (l+m)*(l-m)/((2*l+1)*(2*l-1)))^0.5;
end