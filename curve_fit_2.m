clear,  close, clc
T0=0; h=1.; Tn=5000; [Rung, Newt]=Bfunc;
Ar=1000; a=3;CI=.01;CLS=-3.1;
% -------------------------------------------------------------------------
% Av0R=1/3*[1.,1.e-4,1.e-4,1.1,.1]'; Av0N=[.35,.0001,.0001,.55,.1]';
% E=.1;G=10*E; dU=zeros(3,3,2); dU(6)=G; dU(9+[1 6 9])=[E G -E];
% id=[2 1 1 7.0] ;nid=length(id); kp={{.1}, {1.}, {.1},{.9, .0}};
% for i=1:2
%     for j=1:nid
%         k=j+(i-1)*nid;
%         var={Ar, a,@(t) dU(:,:,i),@(t) [],id(j),CLS,[kp{j}, {CI}],{2,3,1}};
%         [~,Avn2(j,:,i)]=Newt([],Av0N,var{:});
%         [tn(j,i),Avn1(j,:,i),t{i}(:,j),Av{i}(:,:,j)]=Rung(T0,h,Tn,Av0R,var{:});k
%     end
% end
% 
% clr={'k', 'g','r','c'}; nms={'RSC', "FT ", 'SRF','RPR'};
% lst={'--', '-', '-.'}; pk=[1 4 5];
% for i=1:2
%     f=figure(i); clf; f.Color='w'; grid on; hold on
%     f.Position=[705   280   545   435];
%     for j=1:nid
%         for k=1:3
%             Axx="$\rm "+nms{j}+"-a_{"+string(5*(pk(k)-k)+10+k)+"}$";
%             p=plot(G*t{i}(:,j),Av{i}(:,pk(k),j),[clr{j},lst{k}],'DisplayName',Axx,...
%                 'LineWidth',.5);
%         end
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',14);
%     ylabel("$\rm a_{ij}$",'Interpreter','latex','FontSize',14);
%     xlim([0 400]);ylim([0,1]);
%     set(gca,'TickDir','both','GridLineStyle','--','Box','on'); 
%     ax=gca; obj=ax.findobj(); obj2=obj([13,10,7,4,12,9,6,3,11,8,5,2]);
%     legend(obj2,'Location','southoutside','Orientation','horizontal',...
%         'NumColumns',4,'FontSize',10,'Box','off','Interpreter','latex');
% end
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100; err(isnan(err))=0;
% T=array2table(abs([err(:,[1 4 5],1) err(:,[1 4 5],2)]));
% T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.VariableNames=...
%     {'A_11_L1','A_22_L1','A_13_L1','A_11_L2','A_22_L2','A_13_L2'};
% T.Properties.RowNames={'RSC', 'FT', 'SRF','RPR'};
% -------------------------------------------------------------------------
Av0N=[.35,.0001,.0001,.55,.1]'; Av0R=1/3*[1.,1.e-4,1.e-4,1.,.01]';
Ar=1000; a=3;CLS=-4.1; G=1;
dU=zeros(3); dU(6)=G;   
CI={.0165;.0630;.0060}; ck={.988; .965; .900}; CM={.999; 1.010; .900}; 
kpa={1/30; 1/30; 1/20}; ak=num2cell(1-[kpa{:}]'); bk={0.0;0.0;0.0};
bta={[3.842 -17.86  525.  .1168  -5. ]*1e-4;
     [37.28 -169.5  1750. -33.67 -100]*1e-4;
     [4.643 -6.169  190.    9.65  7  ]*1e-4};
for j=1:3, Dk{j,1}=diag([1, ck{j}, 1-ck{j}]); end
val={[ak, bk, CI, CM], [ak, bk, CI, Dk], [kpa, bta]};
id=[7.4, 7.3, 6.1] ;nid=length(id);
for i=1:3
    for j=1:nid
        k=j+nid*(i-1);
        var={Ar, a,@(t) dU,@(t) [],id(j), CLS,val{j}(i,:),{2,4}};
        [~, Avn2(k,:)]=Newt([],Av0N,var{:});[i,j]
        [tn(k),Avn1(k,:),t{i,j},Av{i,j}]=Rung(T0,h,Tn,Av0R,var{:});
    end
end
%
Lns={'-', '--',':'}; pk=[1 4 5];
nms={'iARD-RPR','pARD-RPR', 'iARD-RSC'}; clr={'r','g','c'};
ttl={'40%wt.glass-fiber/PP', '31%wt.carbon-fiber/PP', '40%wt.glass-fiber/nylon'};
for i=1:3
    f=figure(i); clf; f.Color='w'; f.Position=[680,370,640,505];
    grid on; hold on
    % title(['\rm\it' ttl{i}]);
    for j=1:nid
        for k=1:3
            Axx="$\rm "+nms{j}+"-a_{"+string(5*(pk(k)-k)+10+k)+"}$";
            plot(G*t{i,j},Av{i,j}(:,pk(k)),Lns{k},'DisplayName',Axx,...
                'LineWidth',1.,'Color',clr{j});
        end
    end
    xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',14);
    ylabel("$\rm a_{ij}$",'Interpreter','latex','FontSize',14);
    set(gca,'TickDir','both','GridLineStyle','--','Box','on');
    xlim([0 G*Tn]); ylim([0,1]); 
    legend('Location','southoutside','Orientation','horizontal','NumColumns',3,...
        'FontSize',10,'Box','off','Interpreter','latex'); 
end
err=round(Avn1-Avn2,6)./round(Avn1,6)*100; err(isnan(err))=0;
T=array2table(abs(err(:,[1 4 5])));T=varfun(@(x) num2str(x, '%.4f'),T);
T.Properties.VariableNames={'A_11','A_22','A_13'};
Tn=['iARD-RPR';'pARD-RPR'; 'iARD-RSC'];nTn=size(Tn,1);
Tn=[repmat(Tn,3,1),repelem(num2str((1:3)'),nTn,1)];nTn=size(Tn);
T.Properties.RowNames=mat2cell(Tn,ones(nTn(1),1),nTn(2));
%--------------------------------------------------------------------------
