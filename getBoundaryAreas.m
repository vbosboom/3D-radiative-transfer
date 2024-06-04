function Delta = getBoundaryAreas(p,e)
p1 = p(:,e(1,:));
p2 = p(:,e(2,:));
p3 = p(:,e(3,:));

d1 = p2-p1;
d2 = p3-p1;

Delta = sqrt(sum(cross(d1,d2).^2,1));
end