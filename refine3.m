function [pn,en,tn,prol,pc]=refine3(p,e,t)

  [ue,te,ee]=makeedgelist3(t,e);
  pn = 0.5*(p(:,ue(1,:))+p(:,ue(2,:)));
  pn  = [p,pn];

  np = size(p,2);
  te = te+np;
  tt = [t;te];
   
  tn = [tt(5,:), tt(5,:), tt(6,:),  tt(8,:),  tt(9,:), tt(9,:), tt(8,:),  tt(9,:);...
        tt(8,:), tt(6,:), tt(7,:),  tt(9,:),  tt(8,:), tt(7,:), tt(9,:),  tt(6,:);...
        tt(7,:), tt(9,:), tt(10,:), tt(10,:), tt(7,:), tt(6,:), tt(7,:),  tt(7,:);...
        tt(1,:), tt(2,:), tt(3,:),  tt(4,:),  tt(5,:), tt(5,:), tt(10,:), tt(10,:)];
  
  ee = ee+np;
  ee = [e;ee];
  en = [ ee(1,:), ee(2,:), ee(3,:), ee(4,:); ...
         ee(4,:), ee(5,:), ee(6,:), ee(5,:); ...
         ee(6,:), ee(4,:), ee(5,:), ee(6,:)];

%prolongation operator
npn = size(pn,2)-np; %number of new points
kk=(1:npn)';
ue=ue';
ind=[kk,ue(:,1),0.5*ones(npn,1); ...
     kk,ue(:,2),0.5*ones(npn,1)];
prol2=sparse(ind(:,1),ind(:,2),ind(:,3),npn,np);
prol=[speye(np);prol2];

%prolongation for pcw. constants
% col=repmat((1:size(t,2))',1,8);
% row=col + repmat((0:7)*size(t,2),size(t,2),1);
% max(row)
% max(col)
% size(tn,2)
% size(row)
% size(col)
pc=kron(ones(8,1),speye(size(t,2)));
% % sparse(row,col, ones(size(col)),size(tn,2),size(t,2));
% size(pc)
% spy(pc)
% full(pc(1:10,1:10))

      
function [ued,te,ee]=makeedgelist3(t,e)
  
  %% make list of all tet edges 
  ed=sort([t(1,:),t(2,:),t(3,:),t(1,:),t(2,:),t(3,:);...
           t(2,:),t(3,:),t(1,:),t(4,:),t(4,:),t(4,:)]);
  [ued,ie,je]=unique(ed','rows'); ued=ued'; %% ued=ed(:,ie), sed=ued(:,je) 

  %% tet to edge list
  ned=size(ed,2);
  te=reshape(je,ned/6,6)';

  %% make list of all boundary edges
  ed=sort([e(1,:),e(2,:),e(3,:);...
           e(2,:),e(3,:),e(1,:)]);
  [bb,ib,jb]=unique([ed';ued'],'rows');

  %% trig to edge list
  jb=jb(1:length(jb)-length(ie));
  ned=size(ed,2);
  ee=reshape(jb,ned/3,3)';  

