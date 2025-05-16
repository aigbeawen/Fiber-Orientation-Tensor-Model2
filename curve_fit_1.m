clear,  close, clc
[Rung, Newt]=Bfunc;
Av0N=[.3,.0001,.0001,.6,.1]'; Av0R=1/3*[1.,1.e-4,1.e-4,1.,1.e-4]'; 
Ar=1000; a=3; t0=0; ti=0.005; tn=500;
mkr={'o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'};
% -------------------------------------------------------------------------
% id2=[1.0,3.1:.1:3.3, 3.5,3.6, 4.3] ; nid2=length(id2); 
% lgn={'FT','PT','iARD','pARD','WPT','Dz','MRD'};
% clr2=rand(nid2,3); cls=-4.1; dU=@(t) DU(t,1,1);
% for j=1:nid2
%     val=mdl(id2(j));
%     var={Ar, a, dU, @(t) [],id2(j),cls,val,{2,4}};
%            [~,Avn2(j,:)]=Newt([],Av0N,var{:});
%     [tn(j),Avn1(j,:),t(:,j),Av(:,:,j)]=Rung(t0,ti,tn,Av0R,var{:});
% end
% %
% for m=1:2
%     f=figure(m);clf,f.Color='w'; hold on
%     Axx=['\rm a_{' num2str(11*(3-m)) '}'];
%     % title([Axx '-Component'],'FontSize',12);
%     for j=1:nid2
%         mkrj=mkr(mod(j-1,15)+1);
%         % plot(t(:,j),Av(:,3*m-2,j),'Marker',mkrj,'Color',clr2(j,:),...
%         %     'MarkerSize',.5,'LineStyle','-','LineWidth',1);
%         plot(t(:,j),Av(:,3*m-2,j),'Color',clr2(j,:),...
%             'LineStyle','-','LineWidth',.5);
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex');
%     ylabel("$"+Axx+"$",'Interpreter','latex','FontSize',14);
%     xlim([0 100]); legend(lgn{:},'Box','off');
%     set(gca,'GridLineStyle','--','TickDir','both','Box','on');
%     grid on
% end
% %
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100        ;err(isnan(err))=0;
% T=array2table(abs(err(:,[1 4 5],1)));T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.RowNames=lgn;T.Properties.VariableNames={'A_11','A_22','A_13'};
% -------------------------------------------------------------------------
% cls{1}=[1.1,1.2,2.1,2.2,2.3,2.4,2.5,2.6];
% dsn{1}=["HYB_1","HYB_2","ISO","LIN","QDR","SF2","HL1","HL2"];
% 
% cls{2}=[3.0,-1.1,-2.1,-2.3,-2.4,-3.1,-3.2,-3.3,-3.4,-4.1,-4.2,-4.3]; 
% dsn{2}=["IBOF","ORS","ORT","NAT_1","ORW","NAT_2",...
%         "WTZ","LAR32","ORW3","VST","FFLAR4","LAR4"];
% G=1; dU=@(t) DU(t,G,1);ncls(1)=length(cls{1});ncls(2)=length(cls{2});
% %
% for i=1:2
% nclsi=ncls(i);clr2=rand(nclsi,3); 
% for j=1:nclsi
%     var={Ar, a, dU, @(t) [], 1.,cls{i}(j),{1., .01},{2,4}};
%     [~      ,Avn2(j,:,i)                     ]=Newt([]      ,Av0N,var{:});
%     [tn(j,i),Avn1(j,:,i),t(:,j,i),Av(:,:,j,i)]=Rung(t0,ti,tn,Av0R,var{:});
% end
% end
% %
% for i=1:2
%     for j=1:2
%         f=figure(j+2*(i-1)); clf,f.Color='w'; hold on
%         Axx=['\rm a_{' num2str(10+j) '}'];
%         % title([Axx '-Component']);
%         for k=1:ncls
%             mkrk=mkr(mod(k-1,15)+1);
%             plot(G*t(:,k,i),Av(:,j+3,k,i),'Marker',mkrk,'Color',clr2(k,:),...
%                 'MarkerSize',.5,'LineWidth',.5);
%         end
%         xlabel('\it\.{$\gamma$}t','Interpreter','latex');
%         ylabel("$"+Axx+"$",'Interpreter','latex','FontSize',14);
%         xlim([0 500]);  ylim([0,1]-[.2,.7]*(j-1));
%         grid on; set(gca,'Box','on'); set(gca,'Tickdir','both');
%         legend(dsn{i}(:),'Location','southoutside','NumColumns',5,...
%             'Orientation','horizontal','FontSize',10,'Box','off');
%     end
% end
% 
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100; 
% err(isnan(err))=0;err=err(:,[1 4 5],:);
% for j=1:2
%     Tj=array2table(abs(err(1:ncls(j),:,j)));
%     T{j}=varfun(@(x) num2str(x, '%.4f'),Tj);
%     T{j}.Properties.VariableNames={'A_11','A_22','A_12'};
%     T{j}.Properties.RowNames=dsn{j};
% end
%% ------------------------------------------------------------------------
function f=mdl(flg)
flg1=floor(flg); flg2=round(10*(flg-flg1));
switch flg1
    case 1 % FT - SRF
        f={ 1 .0311};
    case 2 % RSC
        f={.2 .0311};
end
switch flg1
case {3 4 6 7}
switch flg2
    case 0 % IRD
        f={.0311}; % 'CI'
    case 1 % ARD
        f={[1.924,58.390,400.,.1168,0.]*1e-4}; % 'bta 1-5'
    case {2 4} % iARD
        f={.0562 .9977};     % ['CI' 'CM'] ;
    case 3 % pARD
        switch flg1
        case 4 % mARD
           %D=zeros(3); D([1 5 9])=[1.,.8,.15];f={.04796 D}; % ['CI' 'D' ] ;
            D=zeros(3); D([1 5 9])=[1.,.7946,.012];f={.0198 D}; % ['CI' 'D' ] ;
        otherwise % pARD
            c=.9868; D=[1,c,1-c]; D=diag(D,0);
            f={.0169 D};      % ['CI' 'D' ] ;
        end
    case 5 % WPT
        f={.0504 .9950}; % ['CI', 'w']
    case 6 % Dz
        n=[0,0,1];
        f={.0258, .0051, n}; % ['CI', 'Dz', 'n']
end
case 5 % NEM
    f={.0311; 1};  % ['CI, U0']
end
switch flg1
    case 6 % (p)ARD-RSC
        f=[{1/30} f]  ; %['kpa' ARD-vars]
    case 7 % X-RPR
        f=[{1-1/30, .01} f]; %['lfa' 'bta', X-vars]
end
end
%% ------------------------------------------------------------------------
function f=DU(t,flg,varargin)
f=zeros(3);
switch flg
    case 1 % Simple shear
        G=varargin{1};  f(6)=G;
    case 2 % Shearing/stretching
        [E,G]=varargin{1:2}; 
        f([1,5,6,9])=[2*E,-E,G,-E];
    case 3 % Uniaxial
        E=varargin{1}; 
        f([1,5,9])=[-E,-E,2*E];
    case 4 % Biaxial
        E=varargin{1};
        f([1,5,9])=[-2*E, E, E];
    case 5 % Shearing/planar stretching
        [E,G]=varargin{1:2};
        f([3,5,9])=[G,E,-E];

end
end