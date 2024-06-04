%We sort the spherical harmonics by the list (l,m) = (0,0), (1,-1) (1,0)
%(1,1) (2,-2), .... This function returns the index of the harmonic (l,m) in
%this list or an empty vector if it does not exist
function indx = momentToDof(l,m)
    indx = [];
    if l>=0 && m<=l && m>=-l %check if l and m are valid number
        if mod(l,2) == 1 %odd moments
            indx = sum(2*(1:2:(l-2))+1)+l+m+1; %conversion formula for index
        end
        if mod(l,2) == 0 %even moments
            indx = sum(2*(0:2:(l-2))+1)+l+m+1; %conversion formula for index
        end
    end
end