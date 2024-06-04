%assemble the weighted mass matrices \int_s f(s) H_l^m(s) \bar(H_k^n(s)) ds
%for A1: f(s) = (sx-i*sy)/2, for A2: f(s) = (sx+i*sy)/2, for A3: f(s) = sz
%for information on coefficients see https://publications.rwth-aachen.de/record/82671/

function [A1,A2,A3] = getAngularA(N)


MomentToDof=@(l,m) momentToDofFull(l,m);


Ne = sum(2*(0:2:N)+1);
No = sum(2*(1:2:N)+1);
A1 = zeros(No,Ne); A2 = zeros(No,Ne); A3 = zeros(No,Ne);

for l=0:2:N
    for m=-l:l
        [k,n,r] = getFirstA1(l,m);  A1(MomentToDof(k,n),MomentToDof(l,m)) = r;
        [k,n,r] = getSecondA1(l,m); A1(MomentToDof(k,n),MomentToDof(l,m)) = r;
        
        [k,n,r] = getFirstA2(l,m);  A2(MomentToDof(k,n),MomentToDof(l,m)) = r;
        [k,n,r] = getSecondA2(l,m); A2(MomentToDof(k,n),MomentToDof(l,m)) = r;
        
        [k,n,r] = getFirstA3(l,m);  A3(MomentToDof(k,n),MomentToDof(l,m)) = r;
        [k,n,r] = getSecondA3(l,m); A3(MomentToDof(k,n),MomentToDof(l,m)) = r;
    end
end
A1 = sparse(A1); A2=sparse(A2); A3=sparse(A3);
end



function [k,n,r] = getFirstA1(l,m)
if l == 0  || m <= -l+1
    k=[]; n =[]; r= [];      
else
    k = l-1;
    n = m-1;
    r = a(l,m);
end

end
function [k,n,r] = getSecondA1(l,m)
k = l+1;
n = m-1;
r = -a(l+1,1-m);
end

function [k,n,r] = getFirstA2(l,m)
if l == 0  || m >= l-1
    k=[]; n =[]; r= [];      
else
    k = l-1;
    n = m+1;
    r = -a(l,-m);
end
end

function [k,n,r] = getSecondA2(l,m)
k = l+1;
n = m+1;
r = a(l+1,m+1);
end

function [k,n,r] = getFirstA3(l,m)
if l == 0  || abs(m) == l
    k=[]; n =[]; r= [];      
else
    k = l-1;
    n = m;
    r = b(l,m);
end
end
function [k,n,r] = getSecondA3(l,m)
k = l+1;
n = m;
r = b(l+1,m);
end

function i = momentToDofFull(l,m)
    if isempty(l) || isempty(m)
        i = [];
    elseif mod(l,2) == 1 % odd moment
        no = sum(2*(1:2:(l-2))+1);
        i = no + l+m+1;
    elseif mod(l,2) == 0 % even moment
        ne = sum(2*(0:2:(l-2))+1);
        i = ne + l+m+1;
    end
end

function erg=a(l,m)
erg=0.5*((l+m)*(l+m-1)/((2*l+1)*(2*l-1)))^0.5;
end

function erg=b(l,m)
erg= ( (l+m)*(l-m)/((2*l+1)*(2*l-1)))^0.5;
end