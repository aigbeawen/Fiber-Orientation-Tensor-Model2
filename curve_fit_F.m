clear,  close, clc
[Rung, Newt]=Bfunc;
Av0R=1/3*[1.,1.e-4,1.e-4,1.,1.e-4]'; Ar=1000; a=3; t0=0; ti=0.01; tn=100;
mkr={'o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'};
lgd={'$SS$','$SUA_1$','$SUA_2$','$UA$','$BA$','$PST_1$','$PST_2$',...
     '$SBA_1$','$SBA_2$','$TA$','$STA_1$','$STA_2$'}; 
% -------------------------------------------------------------------------
id2=1.0; kp=1.; CI=0.01; val2={kp CI}; cls=-2.3; % -2.3, -3.1, -3.2, 3.0
% -------------------------------------------------------------------------
% val1={1;[0.1, 1;10 1]; 1; 1;[0.1, 1;10, 1];[0.2, 1;.5, 1];1;[0.2, 1;.5, 1]};
% nf=length(val1); jj=0;
% for j=1:nf
%     valj=val1{j}; ns(j)=size(valj,1);
%     for k=1:ns(j)
%         jj=jj+1; [dU,Av0N]=DU(j,valj(k,:));
%         var={Ar, a,@(t) dU , @(t) [],id2,cls,val2,{2,4}};
%         [~     , Avn2(jj,:)                   ]=Newt([],Av0N,var{:}); [j,k]
%         [tn(jj), Avn1(jj,:),t(:,jj),Av(:,:,jj)]=Rung(t0,ti,tn,Av0R,var{:});
%     end
% end
% lsty={'-','--','-.'}; clr={'r','g','g','b','c','b','m','m'};
% % 
% for m=1:3
%     f=figure(m);clf;f.Color='w'; grid on; hold on
%     n=-(m^2-6*m+4); jj=0;
%     Axx=['$\rm a_{' num2str(.5*(13*m^2-61*m+92)) '}$'];
%     % title([Axx '-Component']);   
%     for j=1:nf
%         nf2=ns(j);
%         for k=1:nf2
%             jj=jj+1; mkrj=mkr(mod(jj-1,15)+1);
%             % plot(t(:,jj),Av(:,n,jj),'Marker',mkrj,'Color',clr2(jj,:),...
%             %     'MarkerSize',1,'LineWidth',.5,'DisplayName',...
%             %     ['\it' lgd{j} '-' num2str(k)]);
%             plot(t(:,jj),Av(:,n,jj),'Color',clr{j},'LineStyle',lsty{k+(nf2>1)},...
%                 'LineWidth',.5,'DisplayName',['\it' lgd{j} '-' num2str(k)]);
%         end
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',14);
%     ylabel(Axx,'Interpreter','latex','FontSize',14); xlim([0 30]); 
%     legend(lgd{:},'Location','southoutside','Orientation','horizontal',...
%         'Box','off', 'FontSize',10,'Interpreter','latex','NumColumns',4);
%     ax=gca; obj=ax.findobj(); obj2=obj([1,10,9,4,13,12,6,3,8,11,5,2,7]);
%     legend(obj2(2:13)); f.Position=[1075,300,470,430];
%     set(gca,'TickDir','both','GridLineStyle','--','Box','on');
% end
% %
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100        ;
% err(isnan(err))=0; T=array2table(abs(err(:,[1 4 5]))) ;
% T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.VariableNames={'A_11','A_22','A_13'};
% T.Properties.RowNames=lgd;
%--------------------------------------------------------------------------
% nv=21; va=linspace(.01,8,nv)'; vb=logspace(-1.5,2,nv)';
% v1=ones(nv,1); v2=1./va;  v3=1./vb;
% val1={1;[v2 v1]; 1; 1;[v2 v1];[v3 v1];1;[v3 v1]};nf=length(val1); 
% 
% jj=0;
% for j=1:nf
%     valj=val1{j}; ns(j)=size(valj,1);
%     for k=1:ns(j)
%         jj=jj+1; [dU,Av0N]=DU(j,valj(k,:));
%         var={Ar, a,@(t) dU , @(t) [],id2,cls,val2,{2,4}};
%         [~     , Avn2(jj,:)                   ]=Newt([],Av0N,var{:});
%     end
% end
% f=figure(4); f.Color='w'; 
% for j=1:size(Avn2,1)
%     v=sort(eig(v2M(Avn2(j,:),1)),'descend');
%     plot(v(1),v(2),'+'); hold on
% end
% X=[
%     0    1    1    0
%     0  0.5    1  0.5
%     0  0.5    1    0
%     0    1  0.5    0
%     0  0.5    0  0.5];
% for j=1:5
%     line(X(j,1:2),X(j,3:4),'Color','k')
% end
% text(1/3,1/3,'(1/3,1/3)','FontSize',12);
% text(1/2,1/2,'(1/2,1/2)','FontSize',12);
% text(1  ,0  ,'(1,0)'    ,'FontSize',12);
% text(0  ,1  ,'(0,1)'    ,'FontSize',12);
% set(gca,'Box','on','FontSize',14);
% ax=gca; ax.XAxis.TickLabels='';ax.YAxis.TickLabels='';
% xlabel('\it \lambda_{2}');ylabel('\it \lambda_{1}');
% -------------------------------------------------------------------------
function [f, Av0]=DU(flg,val)
f=zeros(3);
switch flg
    case 1 % Simple shear, SS
        G=val(1);  f(6)=G;
        Av0=[.35,1.e-4,1.e-4,.55,1.e-4]';
    case 2 % Shearing/stretching, SUA
        E=val(1); G=val(2); f([1,5,6,9])=[2*E, E, G, -E];
        Av0=[.7,1.e-4,1.e-4,.2,1.e-4]';
    case 3 % Uniaxial, UA
        E=val(1); f([1,5,9])=[-E,-E,2*E];
        Av0=[.1,1.e-4,1.e-4,.1,1.e-4]';
    case 4 % Biaxial, BA
        E=val(1); f([1,5,9])=[E, E,-2*E];
        Av0=[.4,1.e-4,1.e-4,.4,1.e-4]';
    case 5 % Shearing/planar stretching, PST
        E=val(1);G=val(2); f([1, 6, 9])=[E,  G, -E];
        Av0=[.7,1.e-4,1.e-4,.2,.01]';
    case 6 % Balanced shear/biaxial elongation flow, SBA
        E=val(1);G=val(2); f([1,5, 6, 9])=[-2*E, E, G, E];
        Av0=[.1,1.e-4,1.e-4,.8,1.e-4]';
    case 7 % Triaxial, TA
        E=val(1); f([1,5,9])=[E, E, E];
        Av0=[.4,1.e-4,1.e-4,.4,1.e-4]';
    case 8 % Balanced shear/triaxial elongation flow, STA
        E=val(1);G=val(2); f([1,5, 6, 9])=[E, E, G, E];
        Av0=[.4,1.e-4,1.e-4,.4,1.e-4]';
end
end
%%
function A=v2M(Av,flg)
for m=1:2
    for n=m:3
        k=2*(m-1)+n;
        switch flg
            case 1
                A(m,n)=Av(k); A(n,m)=A(m,n);
            case 2
                A(k,1)=Av(m,n);
        end
    end
end
if flg==1, A(3,3)=1-A(1,1)-A(2,2); end
end