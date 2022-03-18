% k_elemu
function [KE]=K_elemu(prenode,postnode,H,Ki,Kgb,k, Kele, ndg)
% PURPOSE : This is a subprogram as for the Stiffness matrix
%
h=postnode-prenode;  
Kd=double(subs(Kele)); 

if ndg==3
%Only for HM Matrix
%Organizing matrix
KE=zeros(6,6);
KE(1:2,:)=Kd(1:2,:);
KE(3,:)=Kd(5,:);
KE(4:5,:)=Kd(3:4,:);
KE(6,:)=Kd(6,:);
KEd=KE;
KE(:,1:2)=KEd(:,1:2);
KE(:,3)=KEd(:,5);
KE(:,4:5)=KEd(:,3:4);
KE(:,6)=KEd(:,6);
else
    KE=Kd;
end

end