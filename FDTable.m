clear, clc
format shortg
[Rung, ~,ddA,dMdA]=BfuncC; v=1; a=3; re=1000; 
% A=rand(3); A=1/2*(A+A'); Ad=diag(A,0); Adev=.1*(A-diag(Ad,0));
% Aext=diag(Ad+1/3*(1-trace(A)),0);A=Aext+Adev;Av=[A(1,1:3), A(2,2:3)];
Av=[0.0622,0.0765,0.0398,0.5521,0.0186];
G=1; E=1; dU=zeros(3); dU([1,5, 6, 9])=[-2*E, E, G, E];
% -------------------------------------------------------------------------
mdl_id=[1.0,3.1:.1:3.3, 3.5,3.6, 5.0, 6.3, 7.2] ; 
mdl_nm={'FT','PT','iARD','pARD','WPT','Dz','NEM', 'pARD-RSC', 'iARD-RPR'};
cls_id=[1.1, 1.2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 3.0,...
       -1.1,-2.1,-2.3,-2.4,-3.1,-3.2,-3.3,-3.4,-4.1,-4.2,-4.3];
cls_nm=["HYB_1","HYB_2","ISO","LIN","QDR","SF2","HL1","HL2",...
        "IBOF","ORS","ORT","NAT_1","ORW","NAT_2",...
        "WTZ","LAR32","ORW3","VST","FFLAR4","LAR4"];
n1=length(mdl_id); n2=length(cls_id); err=zeros(n1,n2);
%
for j1=1:n1
    jmdl=mdl_id(j1);
    for j2=1:n2
        invars={v, a, @(t) dU ,@(t) [], jmdl,cls_id(j2),mdl(jmdl),{2,4}};
        [~, J1]=ddA([],Av, invars{:});
        [~ ,J2]=dMdA(@(x) ddA([], x, invars{:}), Av, 1.e-6, 6);
        J2=reshape(J2,5,5); err(j1,j2)=norm(J1-J2);
    end
end
arr={1:8,9:14, 15:n2};
for j=1:3
Tj=array2table(err(:,arr{j})); Tj=varfun(@(x) num2str(x, '%.4e'),Tj);
Tj.Properties.VariableNames=cls_nm(arr{j});Tj.Properties.RowNames=mdl_nm;
T{j}=Tj;
end
% -------------------------------------------------------------------------
function f=mdl(flg)
flg1=floor(flg); flg2=round(10*(flg-flg1));
switch flg1
    case 1 % FT - SRF
        f={1. .0311}; % ['kp' 'CI']
    case 2 % RSC
        f={.2 .01}; % ['kp' 'CI']
end
switch flg1
    case 3
        switch flg2
            case 0 % IRD
                f={.0311}; % 'CI'
            case 1 % ARD
                f={[1.924,58.390,400.,.1168,0.]*1e-4}; % 'bta 1-5'
            case {2 4} % iARD
                f={.0562 .9977};     % ['CI' 'CM'] ;
            case 3 % pARD
                c=.9868; D=[1,c,1-c]; D=diag(D,0);
                f={.0169 D};      % ['CI' 'D' ] ;
            case 5 % WPT
                f={.0504 .9950}; % ['CI', 'w']
            case 6 % Dz
                n=[0,0,1];
                f={.0258, .0051, n}; % ['CI', 'Dz', 'n']
        end
    case 5 % NEM
        f={.01; .063};  % ['CI, U0']
    case 6 % pARD-RSC
        CI=.027; c=.95; kp=.8; D=[1,c,1-c]; D=diag(D,0);f={CI D};f=[{kp} f];
    case 7 % iARD-RPR
        CI=.063; CM=.995; lfa=.2; bta=.01;
        f={lfa bta, CI, CM}; 
end
end