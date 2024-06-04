% inverse matrix of scattering term applied to odd part of the equation:
% C2i = (kron(Mo,mua*I+mus*(1-theta_l*I))^(-1)
function C2i = getC2Inverse(N,mua,mus,thetaOdd,volTet)

No = sum(2*(1:2:N)+1);
C2i = zeros(length(volTet),No);

indx=1;
for l=1:2:N
    for m=-l:l
        d = volTet.*( mua + (1-thetaOdd(indx)) *mus);
        C2i(:,indx) = 1./d;
        indx=indx+1;
    end
end

end