clear,  close, clc, format shortg
[Rung, Newt]=BfuncR;
v=1.; z=.5; a=3; cls=3.0;  r0=1.; ri=1.e-3;
Av0R=1/3*[1.,0.,0.,1.,0.]'; rn=1000; Av0N=[.2,1.e-5,1.e-5,.6,-.05]';
%
mdl_id=[6.1,7.4, 7.3]; nid=length(mdl_id);
mdl_pr={
        {1/30 [3.842 -17.86  525.  .1168  -5. ]*1e-4} ;
        {1-1/30 0. .0165 .999                       } ; % lfa bta CI CM-iARD-RPR 
        {1-1/30 0. .0165 diag([1. .988 0.012],0)    }}; % lfa bta CI  D-pARD-RPR           
for j=1:nid
    var={v, a, @DU, z, mdl_id(j), cls,mdl_pr{j},{2,4}}   ;
    [rN(j) ,Avn2(j,:)                 ]=Newt(      rn,Av0N,var{:});
    [rR(j) ,Avn1(j,:),r(:,j),Av(:,:,j)]=Rung(r0,ri,rn,Av0R,var{:});
end
%
lgd={'ARD-RSC','iARD-RPR','pARD-RPR'};
clr={'r','g','k'}; mkr={'+','x','p'};lns={'--',':','-.'};
clf; f=figure(1); f.Color='w';  hold on
for j=1:nid
    plot(r(:,j),Av(:,[1 4 5],j),'Color',clr{j},'Marker',mkr{j},...
        'MarkerSize',.5,'LineWidth',.5,'LineStyle',lns{j});
end
grid on; title('\rm\it z/h=0.5')
legend('ARD-RSC','','','iARD-RPR','','','pARD-RPR','Location','best'); 
xlabel('\it r/h');ylabel('\it A_{xx}');
set(gca,"Box","on");set(gca,"TickDir",'both');xlim([0 rn]);ylim([-.2 .7]);
%
xj=[.6 .5;.5 .6;.5 .6]; yj=[.78 .87;.73 .61;.35 .27];
for j=1:3
    xx=num2str(11*j-10*(j-1)*(j-2));
    annotation("textarrow",xj(j,:),yj(j,:),'String',['\it A_{' xx '}']);
end
err=abs((Avn2-Avn1)./Avn1)*100; err=err(:,[1 4 5]);T=array2table(err);
T=varfun(@(x) num2str(x, '%.4f'),T);
T.Properties.VariableNames={'A_11','A_22','A_13'};
T.Properties.RowNames=lgd;
%--------------------------------------------------------------------------
function [f,df]=DU(r,z)
Q=100; b=1.5e-1; sz=sign(z);
vb=(3*Q)/(8*pi*b^2); vr=1/r*(1-z^2) ; vz=sz*z/r; 
f =zeros(3);   f([5,1,6])=vb/b*[-vr/r,vr/r, -2*vz];
if nargout>1
df=zeros(3);  df([5,1,6])=(vb/b^2)*(-2/r)*[-vr/r,vr/r, -vz];
end
end