clear,  close, clc
[Rung, Newt]=BfuncR; 
a=3; cls=3.0; Av0=1/3*[1.,0.,0.,1.,0.]'; tol=1.e-2;
v=1.;  r0=1; rn=30; nz=31; z=linspace(-1.+tol,1.-tol,nz);
% 
%
mdl_id=[6.1,7.4, 7.3]; nid=length(mdl_id);
mdl_pr={
        {1/30 [3.842 -17.86  525.  .1168  -5. ]*1e-4} ; % kpa bta(i)   - ARD-RSC
        {1-1/30 0. .0165 .999                       } ; % lfa bta CI CM-iARD-RPR 
        {1-1/30 0. .0165 diag([1. .988 0.012],0)    }}; % lfa bta CI  D-pARD-RPR  
%
for i=1:nid
    for j=1:nz
        var={ v, a, @DU, z(j),mdl_id(i),cls,mdl_pr{i},{2,4}};
        [rRK(j,i), AvRK(j,:,i)]=Rung(r0,.001,rn,Av0,var{:});[i j]
    end
end
%
f=figure(1);clf, f.Color='w'; hold on
Lns={'r--', 'c:','k-.'};nms={'ARD-RSC','iARD-RPR','pARD-RPR'};
for i=1:4
    s=1/2*(i-1)*(i-2); k=s/3*(i-3);
    for j=1:3
        ydat=(1-k)*(s+(2-i)*AvRK(:,1,j)+(i*(1-s)-1)*AvRK(:,4,j))+k*AvRK(:,5,j);
        plot(z,ydat,Lns{j},'DisplayName',['\it' nms{j}]);
    end
end
p=gca; p.XTick=linspace(-1.,1.,9);p.YTick=linspace(-.2, 1,9);
p.XAxisLocation="bottom";p.YAxisLocation="origin";p.Box='on';
grid(p,"on"); xlabel(p,'\it z/h');ylabel(p,'\it A_{xx}');
p.XAxis.TickDirection="out"; p.YAxis.TickDirection="both";
p.YAxis.TickLabelFormat='%.2f';p.FontSize=7;
legend(p,nms{:},'FontSize',8,'Location','northeast','Box','off');
xlim(p,[-1 1]);ylim(p,[-0.2, 1]);
%
aX=[.38 .325;.65 .575;.65 .6;.75 .7];
aY=[0.72 0.63;0.7 0.615; 0.4 0.28; 0.3 0.22];
for j=1:4
    s=1/6*(j-1)*(j-2)*(j-3); Ann=num2str(11*j*(1-s)+13*s);
    annotation('textarrow',aX(j,:),aY(j,:),'String',['\it A_{' Ann '}'],...
        'FontSize',7,'HeadStyle','deltoid','HeadWidth',5);
end
%--------------------------------------------------------------------------
function [f,df]=DU(r,z)
Q=100; b=1.5e-1; sz=sign(z);
vb=(3*Q)/(8*pi*b^2); vr=1/r*(1-z^2) ; vz=sz*z/r; 
f =zeros(3);   f([5,1,6])=vb/b*[-vr/r,vr/r, -2*vz];
if nargout>1
df=zeros(3);  df([5,1,6])=(vb/b^2)*(-2/r)*[-vr/r,vr/r, -vz];
end
end