% m_elemu
function [ME]=m_elemu(prenode,postnode,Lam, Mele, ndg, I)
% PURPOSE : This is a subprogram for  Mass matrix
%
%  M=1/420[ 156  22*lb   54     -13*lb
%                4*lb^2  13*lb  -3*lb^2  
%                        156    -22*lb
%           symmetric            4*lb^2];      
%

    h=postnode-prenode; 
    Md=zeros(ndg*2,ndg*2);
    Md(1:4,1:4)=double(subs(Mele)); 
    
if ndg==3
%Only for HM Matrix
  %Organizing matrix
ME=zeros(6,6);
ME(1:2,:)=Md(1:2,:);
ME(3,:)=Md(5,:);
ME(4:5,:)=Md(3:4,:);
ME(6,:)=Md(6,:);
MEd=ME;
ME(:,1:2)=MEd(:,1:2);
ME(:,3)=MEd(:,5);
ME(:,4:5)=MEd(:,3:4);
ME(:,6)=MEd(:,6);
else
    ME=Md;
end
  end