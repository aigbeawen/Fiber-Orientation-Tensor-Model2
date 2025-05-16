function varargout=BfuncC
varargout={@Rung, @Newt,@ddA,@dMdA};
end
%%
function [tn, Avn,t,Av]=Rung(t0,ti,tn,Av0,re,varargin)
v=(re^2-1)/(re^2+1); ztol=5;
% opts = odeset('Events',@dzero);
% [t, Av]=ode45(@(t,y) ddA(t,y,v,varargin{:}),t0:ti:tn,Av0,opts);
[t, Av]=ode45(@(t,y) ddA(t,y,v, varargin{:}),t0:ti:tn,Av0);
tn=t(end); Avn=Av(end,:); 
    function [pos,ter,dir] = dzero(t,y)
        pos = round(norm(ddA(t,y, v,varargin{:}),2),ztol);
        ter = 1;
        dir = 0;
    end
end
%%
function [t,Av]=Newt(t,Av,re,varargin)
v=(re^2-1)/(re^2+1);
err=1;
while err>1e-5
    [R, J]=ddA(t,Av, v, varargin{:});
%     [R,J2]=dMdA(@(x) ddA(t,x, v, varargin{:}), Av, 1.e-5, 2);
%     J2=reshape(J2,5,5); norm(J-J2);
    Av=Av-J\R; err=norm(R); 
end
%
end
% upper case for matrix, lower case for scalar, v is lambda
%%
function [R, J]=ddA(t,Av, varargin)
[v, a, DU, DUr, flg1,flg2,var1,var2]=varargin{:};
dU=DU(t); dUr=DUr(t);
%
A=v2M(Av,1); I=eye(3);  
if nargout>1
    dA=DA(I); s3=size(dA,3); dI=zeros(3,3,s3);
end
%
if nargout>1
    [F,dF]=mdlR(A, I, dA, dI, v, a, dU, dUr,flg1,flg2,var1,var2);
    R=[F(1,:) F(2,2:3)]';J=permute([dF(1,:,:) dF(2,2:3,:)],[2,3,1]); 
else
    F     =mdlR(A, I, [], [], v, a, dU,  [],flg1,flg2,var1,var2);
    R=[F(1,:) F(2,2:3)]';
end
if nargout>1
    % Rank Deficient - Normalization
    nc=size(J,2); nr=rank(J);
    if (nc-nr)>0
        if nargout>1
            Rn=norm(R);  Jn=1/Rn*(R'*J);  R =[Rn; R]; J= [Jn;J] ;
        end
    end
end
%
end
%%
function [F,dF]=mdlR(A, I, dA, dI, v, a, dU, dUr, flg1,flg2,var1,var2)
%
flg1_1=floor(flg1); flg1_2=round(10*(flg1-flg1_1));
bol1=(nargout>1); bol2=isempty(dUr); bol=logical(abs(bol1-bol2));
%
switch flg1_1
    case {1 2 6} % SRF RSC
        kp=var1{1};
        switch flg1_1
            case {1 2} % IRD
                CI=var1{2}; var1=[];
            case 6 % ARD
                var1=var1(2:end);
        end
    case 5 % NEM
        [U0, CI]=var1{1:2}; var1=[];
end
switch flg1_1
    case 7 % RPR
        [ak, bk]=var1{1:2}; 
        switch flg1_2
            case 0
                CI=var1{3}; var1=[];
            otherwise
                var1=var1(3:end);
        end
            
end
%
W =1/2*(dU'-dU); Y=1/2*(dU+dU');  y=sqrt(2*sum(Y.* Y' ,'all'));
if bol1
    s3 =size(dA,3);
    if ~bol2
    dWr=1/2*(dUr'-dUr); dYr=1/2*(dUr+dUr'); dyr=2/y*sum(Y.*dYr','all');
    end
end
%--------------------------------------------------------------------------
switch flg1_2
    case 4
        Yin=dU;
    otherwise
        Yin=Y; 
end
if bol
    varin_1={A,I,[],[],Yin,[],flg1,flg2,var1,var2,a};
    switch flg1_1
        case {1 7}
            if flg1_2==0
                varin_1{7}=1;  YA4 =BA4(varin_1{:});
            end
    end
    switch flg1_1
        case 2
            [YA4, YLMA4]=BA4(varin_1{:});
        case {3 5 7}
            if flg1_2>0 || flg1_1==5
                [YA4, CA4,C]=BA4(varin_1{:});
            end
        case 4
            [YA4, C    ]=BA4(varin_1{:});
        case 6
            [YA4, YLMA4, CA4, CM4, CLMA4, C]=BA4(varin_1{:});
    end
end
if bol1
    if ~bol2
    switch flg1_2
        case 4
            dYin=dUr;
        otherwise
            dYin=dYr;
    end
    else
        dYin=[];
    end
    varin_2={A,I,dA,dI,Yin,dYin,flg1,flg2,var1,var2,a};
    switch flg1_1
        case {1 7}
            if flg1_2==0
                varin_2{7}=1;  [dYrA4, YdA4 ]=BA4(varin_2{:});
            end
    end
    switch flg1_1
        case 2
            [dYrA4, dYrLMA4,YdA4, YdLMA4]=BA4(varin_2{:});
        case {3 5 7}
            if flg1_2>0 || flg1_1==5
                [dYrA4, dCrA4, C, YdA4,dCA4, dC]=BA4(varin_2{:});
            end
        case 4
            [dYrA4, C, YdA4, dC]=BA4(varin_2{:});
        case 6
            [dYrA4, dYrLMA4, dCrA4, dCrM4, dCrLMA4,  C, ...
              YdA4,  YdLMA4,  dCA4,  dCM4,  dCLMA4, dC]=BA4(varin_2{:});
    end
    if bol2
        YA4=dYrA4;
        switch flg1_1
            case {2 6}
                YLMA4=dYrLMA4;
                switch flg1_1
                    case 6
                        CM4=dCrM4; CLMA4=dCrLMA4;
                end
        end
        switch flg1_1
            case {3 5 6 7}
                if flg1_2>0 || flg1_1==5
                    CA4=dCrA4;
                end
        end
    else
        switch flg1_1
            case {3 4 5 6 7}
                if flg1_2>0 || flg1_1==5
                    dCr=dC(:,:,1); dC=dC(:,:,2:end);
                end
        end
    end
end
% -------------------------------------------------------------------------
Ah=(W*A-A*W) + v*(Y*A+A*Y-2*YA4); 
if bol1
for m=1:s3
    dAh(:,:,m)=(W*dA(:,:,m)-dA(:,:,m)*W)+...
             v*(Y*dA(:,:,m)+dA(:,:,m)*Y)-2*v*YdA4(:,:,m);
end
if ~bol2
    dAhr=(dWr*A-A*dWr)+ v*(dYr*A+A*dYr-2*dYrA4) ;
end
end
%
switch flg1_2
case 0
    Ad=2*y*CI*(I-a*A);
    if bol1
        for m=1:s3
            dAd(:,:,m)=-2*CI*y*a*dA(:,:,m);
        end
        if ~bol2
            dAdr=2*CI*dyr*(I-a*A);
        end
    end
end
% -------------------------------------------------------------------------
switch flg1_1
case 1 % SRF
    FT=Ah+Ad;  F=kp*FT; 
    if bol1
        dFT=dAh+dAd;
        if ~bol2
            dFTr=dAhr+dAdr; dFT=cat(3,dFTr,dFT);
        end
        dF=kp*dFT;
    end
case 2 % RSC
    FT=Ah+Ad;  Arsc=2*v*YLMA4+Ad; F=FT-(1-kp)*Arsc;
    if bol1
        dFT =dAh +dAd ; dArsc =2*v*YdLMA4 +dAd ; dF =dFT -(1-kp)*dArsc ;
        if ~bol2
            dFTr=dAhr+dAdr; dArscr=2*v*dYrLMA4+dAdr; dFTr=dFTr-(1-kp)*dArscr;
            dF=cat(3,dFTr,dF);
        end
    end
case {3 7} % ARD , X-RPR
    if flg1_2>0
    Aard=y*(2*C-2*trace(C)*A-5*(C*A+A*C)+10*CA4);
    if bol1
        if ~bol2
            dAardr=dyr/y*Aard+y*(2*dCr-2*trace(dCr)*A-5*(dCr*A+A*dCr)+10*dCrA4);
        end
        %
        for m=1:s3
            dAard(:,:,m)=y*(2*dC(:,:,m)-2*(trace(dC(:,:,m))*A+trace(C)*dA(:,:,m))+...
                -5*(dC(:,:,m)*A+C*dA(:,:,m)+dA(:,:,m)*C+A*dC(:,:,m))+10*dCA4(:,:,m));
        end
    end
    end
    switch flg1_1
        case 3 % ARD
            F=Ah+Aard; 
            if bol1 
                dF=dAh+dAard; 
                if ~bol2
                    dFr=dAhr+dAardr;  dF=cat(3,dFr,dF);
                end
            end
        case 7 % X-RPR
            nk=[1 2 3];
            switch flg1_2
                case 0    % IRD
                    F=Ah+Ad; 
                    if bol1
                        dF=dAh+dAd; if ~bol2, dFr=dAhr+dAdr; end
                    end
                otherwise % ARD
                    F=Ah+Aard; 
                    if bol1
                        dF=dAh+dAard; if ~bol2, dFr=dAhr+dAardr; end
                    end
            end
            % RPR Correction 
            if bol1 %
                [~, Q, ~, dQ]=deigvdA(A,I,dA,dI,var2{:});
                eigv=diag(Q\F*Q,0);
                for m=1:s3
                    dv(:,m)=diag(-(Q\dQ(:,:,m)/Q)*F*Q+...
                        Q\dF(:,:,m)*Q+Q\F*dQ(:,:,m),0);
                end
            else
                [~, Q]=deigvdA(A,I,[],[],var2{:}); eigv=diag(Q\F*Q,0);
            end
            eigvt=transpose(eigv);
            bol=1-I; V_iok=ak*(eigv-bk*(eigv.^2+2*prod(eigvt.*bol+I,2)));
            V_iok=diag(V_iok,0); A_iok=-Q*V_iok/Q; 
            F=F+A_iok;
            if bol1
                for m=1:3
                    jk=nk(nk~=m); kj=6-m-jk;
                    dV_iok(m,:)=ak*(dv(m,:)-2*bk*(eigv(m).*dv(m,:)+...
                        sum(eigv(jk).*dv(kj,:),1)));
                end
                dV_iok=permute(dV_iok,[1 3 2]).*repmat(I,1,1,s3);
                for m=1:s3
                    dA_iok(:,:,m)=-(dQ(:,:,m)*V_iok/Q+Q*dV_iok(:,:,m)/Q-...
                        Q*V_iok*(Q\dQ(:,:,m)/Q));
                end
                dF=dF+dA_iok; if ~bol2, dF=cat(3,dFr,dF); end
            end
    end
case 4 % MRD
    Amrd=2*y*(C-trace(C)*A); F=Ah+Amrd;
    if bol1
        for m=1:s3
            dAmrd(:,:,m)=2*y*(dC(:,:,m)-trace(dC(:,:,m))*A-trace(C)*dA(:,:,m));
        end
        dF=dAh+dAmrd; 
        if ~bol2 
            dAmrdr=2*(dyr*(C-trace(C)*A)+y*(dCr-trace(dCr)*A));
            dFr=dAhr+dAmrdr; dF=cat(3,dFr,dF);
        end
    end
case 5 % NEM
    AA4=CA4; Anem=y*U0*(A*A-AA4); F=Ah+.5*Ad+Anem;
    if bol1
        dAA4=dCA4;
        for m=1:s3
            dAnem(:,:,m)=y*U0*(dA(:,:,m)*A+A*dA(:,:,m)-dAA4(:,:,m));
        end
        dF=dAh+.5*dAd+dAnem;
        if ~bol2 
            dAnemr=dyr/y*Anem; dFr=dAhr+.5*dAdr+dAnemr; dF=cat(3,dFr,dF);
        end
    end
case 6 % ARD-RSC
    Arsc=-2*v*YLMA4;
    Aard_rsc=y*(2*(C-(1-kp)*CM4)-2*kp*trace(C)*A-5*(C*A+A*C)+...
         10*(CA4+(1-kp)*CLMA4));
    F=Ah+(1-kp)*Arsc+Aard_rsc;
    if bol1
        for m=1:s3
            dArsc(:,:,m)=-2*v*YdLMA4(:,:,m);
            dAard_rsc(:,:,m)=y*(...
                2*(dC(:,:,m)-(1-kp)*dCM4(:,:,m))-...
                2*kp*(trace(dC(:,:,m))*A+trace(C)*dA(:,:,m))-...
                5*(dC(:,:,m)*A+C*dA(:,:,m)+dA(:,:,m)*C+A*dC(:,:,m))+...
                10*(dCA4(:,:,m)+(1-kp)*dCLMA4(:,:,m)));
        end
        dF=dAh+(1-kp)*dArsc+dAard_rsc;
        if ~bol2
            dArscr=-2*v*dYrLMA4;
            dAard_rscr=dyr/y*Aard_rsc+y*(2*(dCr-(1-kp)*dCrM4)-...
                2*kp*trace(dCr)*A-5*(dCr*A+A*dCr)+10*(dCrA4+(1-kp)*dCrLMA4));
            dFr=dAhr+(1-kp)*dArscr+dAard_rscr; dF=cat(3,dFr,dF);
        end
    end
end
%
end
%%
function varargout=BA4(A,I,dA ,dI ,Y, dY,flg1,flg2,var1,var2,a)
%
flg1_1=floor(flg1); flg1_2=round(10*(flg1-flg1_1));
outflg=[1 2 3 2 3 6 3]; nout=outflg(flg1_1);
%
var3=var2;
switch floor(flg2)
    case {0 1}
        var3={a};
end
%
YA4=zeros(3); 
if nargout>nout
    s3=size(dA,3); YdA4=zeros(3,3,s3);
    [A4, dA4 ]=aijkl(A,I,dA,dI,flg2,var3{:});
else
     A4       =aijkl(A,I,[],[],flg2,var3{:});
end
% Initialization
switch flg1_1
    case {2 6}
        if nargout>nout
             YdLMA4=zeros(3,3,s3);
            [L4,M4, dL4, dM4]=LM(A,I,dA,dI,var2{:});
            [MA4,dMA4]=fun_MA(M4,A4, dM4,dA4);
        else
             YLMA4=zeros(3,3); 
            [L4,M4]=LM(A,I,[],[],var2{:}); MA4=fun_MA(M4,A4);
        end
        switch flg1_1
            case 6
                CM4=zeros(3,3); CLMA4=zeros(3,3);
                if nargout>nout
                    dCM4=zeros(3,3,s3); dCLMA4=zeros(3,3,s3);
                end
        end
end
switch flg1_1
    case {3 5 6 7}
        CA4=zeros(3,3); if nargout>nout, dCA4=zeros(3,3,s3); end
end
% 
switch flg1_1
    case {3 4 6 7} % Spatial Diffusion Tensor
        if nargout>nout
            [C, dC]=Cij(A,I,dA,dI,Y,dY,flg1_2,var1,var2);
            if ~isempty(dY)
                switch flg1_2
                    case {1 2 4}
                        dCr=dC(:,:,1); dC=dC(:,:,2:end);
                end
            end
        else
            C          =Cij(A,I,[],[],Y,[],flg1_2,var1,var2);
        end
        switch flg1_2
            case 4
                L=Y; Y=1/2*(L+L');
                if nargout>nout
                    if ~isempty(dY)
                        dL=dY; dY=1/2*(dL+dL');
                    end
                end
        end
    case 5
        C=A; dC=dA; if ~isempty(dY), dCr=zeros(3); end
end
% -------------------------------------------------------------------------
if ~isempty(dY), Yr=dY; else,  Yr=Y ; end
%
switch flg1_1
    case {3 5 6 7}
        if ~isempty(dY), Cr=dCr; else,Cr=C  ;end
end
%
for i=1:3
    for j=1:3
        YA4(i,j)=0  ; if nargout>nout, YdA4(i,j,1:s3)=0;   end
        switch flg1_1
            case {2 6}
                YLMA4(i,j)=0 ;
                if nargout>nout
                    YdLMA4(i,j,1:s3)=0;
                end
                switch flg1_1
                    case 6
                        CLMA4(i,j)=0; CM4(i,j)=0      ;
                        if nargout>nout
                            dCLMA4(i,j,1:s3)=0; dCM4(i,j,1:s3)=0  ;
                        end
                end
        end
        %
        switch flg1_1
            case {3 5 6 7}
                CA4(i,j)=0 ;  if nargout>nout, dCA4(i,j,1:s3)=0  ; end
        end
        %
        for k=1:3
            for l=1:3
                %
                YA4(i,j) = YA4(i,j)+ Yr(k,l)*A4(i,j,k,l);
                if nargout>nout
                    for n=1:s3
                        YdA4(i,j,n) = YdA4(i,j,n)+Y(k,l)*dA4(i,j,k,l,n);
                    end
                end
                %
                switch flg1_1
                    case {2 6}
                        YLMA4(i,j) = YLMA4(i,j)+...
                            Yr(k,l)*(L4(i,j,k,l)-MA4(i,j,k,l));
                        if nargout>nout
                            for n=1:s3
                                YdLMA4(i,j,n) = YdLMA4(i,j,n)+...
                                    Y(k,l)*(dL4(i,j,k,l,n)-dMA4(i,j,k,l,n));
                            end
                        end
                        %
                        switch flg1_1
                            case 6
                                CLMA4(i,j) = CLMA4(i,j)+...
                                    Cr(k,l)*(L4(i,j,k,l)-MA4(i,j,k,l));
                                CM4(i,j) = CM4(i,j)+Cr(k,l)*M4(i,j,k,l);
                                if nargout>nout
                                    for n=1:s3
                                        dCLMA4(i,j,n) = dCLMA4(i,j,n)+...
                                            C(k,l  )*(dL4(i,j,k,l,n)-dMA4(i,j,k,l,n))+...
                                            dC(k,l,n)*( L4(i,j,k,l  )- MA4(i,j,k,l  ));
                                        dCM4(i,j,n) = dCM4(i,j,n)+...
                                            C(k,l)*dM4(i,j,k,l,n)+dC(k,l,n)*M4(i,j,k,l);
                                    end
                                end
                        end
                end
                %
                switch flg1_1
                    case {3 5 6 7}
                        CA4(i,j) = CA4(i,j)+Cr(k,l)*A4(i,j,k,l);
                        if nargout>nout
                            for n=1:s3
                                dCA4(i,j,n) = dCA4(i,j,n)+...
                                    C(k,l)*dA4(i,j,k,l,n)+dC(k,l,n)*A4(i,j,k,l);
                            end
                        end
                end
            end
        end
    end
end
%
switch flg1_1
    case 1  
        varargout={YA4};
        if nargout>nout
            varargout=[varargout {YdA4}];
        end
    case 2 
        varargout={YA4, YLMA4};
        if nargout>nout
            varargout=[varargout {YdA4 YdLMA4}];
        end
    case {3 5 7}
        % AA4=CA4 - 5
        varargout={YA4, CA4, C};
        if nargout>nout
            if ~isempty(dY), dC=cat(3,Cr,dC); end
            varargout=[varargout {YdA4, dCA4, dC}];
        end
    case 4
        varargout={YA4, C};
        if nargout>nout
            if ~isempty(dY), dC=cat(3,Cr,dC); end
            varargout=[varargout {YdA4, dC}];
        end
    case 6
        varargout={YA4, YLMA4, CA4, CM4, CLMA4, C};
        if nargout>nout
            if ~isempty(dY), dC=cat(3,Cr,dC); end
            varargout=[varargout ...
                {YdA4 YdLMA4, dCA4, dCM4, dCLMA4, dC}];
        end
end
end
%% Spatial Diffusion Tensor
function [C, dC]=Cij(A,I,dA,dI,Y,dY,flg,var1,var2)
if nargout>1
    s3=size(dA,3); dC=zeros(3,3,s3); 
end
switch flg
    case {1 2}
        y=sqrt(2*sum(Y.*Y','all')); 
        if nargout>1, if ~isempty(dY), dy=2/y*sum(Y.*dY','all');end; end
        switch flg
            case 1 % ARD
                b=var1{1}; C=b(1)*I+b(2)*A+b(3)*(A*A')+b(4)*Y/y+b(5)*(Y*Y')/y^2;
                if nargout>1
                    dC =b(2)*dA;
                    for m=1:s3
                        dC(:,:,m)=dC(:,:,m)+b(3)*(A*dA(:,:,m)'+dA(:,:,m)*A');
                    end
                    if ~isempty(dY)
                    dCr=b(4)*(y*dY-Y*dy)/y^2+...
                        2*b(5)*(y*(Y*dY'+dY*Y')-2*(Y*Y')*dy)/y^3;
                    end
                end
            case 2 % iARD
                [CI, CM]=var1{1:2};   C=CI*(I-4*CM*(Y*Y')/y^2);
                if nargout>1
                    if ~isempty(dY)
                    dCr=-4*CI*CM*(y*(Y*dY'+dY*Y')-2*(Y*Y')*dy)/y^3;
                    end
                end
        end
    case 3 % pARD
        [CI,D]=var1{1:2};
        if nargout>1
            [~, Q,~, dQ]=deigvdA(A,I,dA,dI,var2{:});
        else
            [~, Q]=deigvdA(A,I,[],[],var2{:});
        end
        C=CI*Q*D/Q;
        if nargout>1
            for m=1:s3
                dC(:,:,m)=CI*(dQ(:,:,m)*D/Q-Q*D*(Q\dQ(:,:,m)/Q));
            end
            if ~isempty(dY), dCr=zeros(3); end
        end
    case 4 % iARD
        [CI,CM]=var1{1:2}; 
        L=Y; lb=sum(L.*L,'all'); Lb=(L*L')/lb; C=CI*(I-CM*Lb);
        if nargout>1
            if ~isempty(dY)
            dL=dY;    dlbr=sum(2*L.*dL,'all');
            dLbr=(lb*(dL*L'+L*dL')-(L*L')*dlbr)/lb^2;
            dCr=-CI*CM*dLbr;
            end
        end
    case 5 % WPT
        [CI, w]=var1{1:2}; C=CI*((1-w)*I+w*(A*A'));
        if nargout>1
            for m=1:s3
                dC(:,:,m)=CI*w*(A*dA(:,:,m)'+dA(:,:,m)*A');
            end
            if ~isempty(dY), dCr=zeros(3); end
        end
    case 6 % Dz
        [CI, Dz, n]=var1{1:3}; C=CI*(I-(1-Dz)*(n'*n));
        if nargout>1
            if ~isempty(dY), dCr=zeros(3); end
        end
end
if nargout>1, if ~isempty(dY), dC=cat(3,dCr,dC); end, end
end
%%
function [MA,dMA]=fun_MA(M,A,dM,dA)
MA=zeros(3,3,3,3); 
if nargout>1, sz=size(dA); s3=sz(end); dMA=zeros(3,3,3,3,s3); end
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                MA(i,j,k,l)=0; 
                if nargout>1, dMA(i,j,k,l,:)=0;end
                for m=1:3
                    for n=1:3
                        MA(i,j,k,l)=MA(i,j,k,l)+M(i,j,m,n)*A(m,n,k,l);
                        if nargout>1
                        dMA(i,j,k,l,:)=dMA(i,j,k,l,:)+...
                            dM(i,j,m,n,:).*A(m,n,k,l)+M(i,j,m,n).*dA(m,n,k,l,:);
                        end
                    end
                end
            end
        end
    end
end
end
%%
function [L,M, dL, dM]=LM(A,I,dA,dI,varargin)
%
M=zeros(3,3,3,3)  ; L=zeros(3,3,3,3);
if nargout>2
    sz=size(dA); s3=sz(end); dM=zeros(3,3,3,3,s3);dL=zeros(3,3,3,3,s3);
    [egv,Q , degv, dQ]=deigvdA(A,I,dA,dI,varargin{:});
else
    [egv,Q]=deigvdA(A,I,[],[],varargin{:});
end
for m=1:3
    Mm =zeros(3,3,3,3)  ; Lm=zeros(3,3,3,3);
    if nargout>2
        dMm=zeros(3,3,3,3,s3);dLm=zeros(3,3,3,3,s3);
    end
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    f=[i j k l];
                    Q4=1; 
                    if nargout>2, dQ4=zeros(1,1,s3); end
                    for r =1:4
                        Q4=Q4*Q(f(r),m);
                        if nargout>2
                            dQr=dQ(f(r),m,:);
                            for s=1:4
                                if r~=s
                                    dQr=dQr.*Q(f(s),m);
                                end
                            end
                            dQ4=dQ4+dQr;
                        end
                    end
                    Mm(i,j,k,l)= Q4; Lm(i,j,k,l)=egv(m)*Q4;
                    if nargout>2
                        dMm(i,j,k,l,:)=dQ4;
                        dLm(i,j,k,l,:)=egv(m)*dQ4+permute(degv(m,:),[1 3 2]).*Q4;
                    end
                end
            end
        end
    end
    M=M + Mm ; L=L + Lm; 
    if nargout>2, dM=dM + dMm ; dL=dL + dLm; end
end
end
%%
function [F2, dF2]=aijkl(A,I,dA,dI,flg,varargin)
flgs=split(num2str(flg),'.');flg1=str2double(flgs{1}); 
if length(flgs)>1, flg2=str2double(flgs{2});end
%
if nargout>1
    s3=size(dA,3);
end
%
switch sign(flg)
    case {0 1}
        a=varargin{1};
    case -1
        if nargout>1
            [eigv,Q,deigv, dQ]=deigvdA(A,I,dA,dI,varargin{:});
        else
            [eigv, Q]=deigvdA(A,I,[],[],varargin{:});
        end

end
%
switch sign(flg)
    case  {0, 1}
        switch flg1
            case {0 1}
                if nargout>1
                    [F2_1, dF2_1]=gencls(A,dA,3); % quad
                    [F2_2, dF2_2]=gencls(A,dA,2); % lin
                else
                    F2_1=gencls(A,[],3); % quad
                    F2_2=gencls(A,[],2); % lin
                end
                switch flg1
                    case 0
                        f=flg; F2=f*F2_1+(1-f)*F2_2;
                        if nargout>1, dF2=f*dF2_1+(1-f)*dF2_2; end
                    case 1
                        switch flg2
                            case 1 % 1.1
                                a1=a/(a-1); b1=1/(a-1);
                                ysq=sum(A.*A','all'); f=a1*ysq-b1;
                                if nargout>1
                                    df=a1*sum(A.*permute(dA,[2 1 3])+dA.*A',[1 2]);
                                end
                            case 2 % 1.2
                                detj=f_detj(A); f=1-a^a*detj;
                                if nargout>1
                                    [~,ddetj]=f_detj(A,dA); df=-a^a*ddetj;
                                end
                        end
                        F2=f*F2_1+(1-f)*F2_2;
                    if nargout>1
                        dF2=f*dF2_1+(1-f)*dF2_2;
                        for m=1:s3
                            dF2(:,:,:,:,m)=dF2(:,:,:,:,m)+df(m)*(F2_1-F2_2);
                        end
                    end
                end
                %
            case 2 % General permutation closure
                if nargout>1
                    [F2, dF2]=gencls(A,dA,flg2);
                else
                    F2=gencls(A,[],flg2);
                end
            case 3 % IBOF
                if nargout>1
                    [F2, dF2]=A4_IBOF(A,I,dA,dI,varargin{:});
                else
                    F2=A4_IBOF(A,I,[],[],varargin{:});
                end
            otherwise %

        end
        %
    case -1
        if nargout>1
            [av,dav]=binom(eigv, -flg1); av=av'; dav=dav';
            dav=dav(:,1)*deigv(1,:)+dav(:,2)*deigv(2,:);
        else
            av=binom(eigv, -flg1); av=av';
        end
        %
        switch flg1
            case -3 % Cubic
                switch flg2
                    case {2 3} % Rational ellipsoid
                        if nargout>1
                            [av2,dav2]=binom(eigv, -flg1+1); av2=av2'; dav2=dav2';
                            dav2=dav2(:,1)*deigv(1,:)+dav2(:,2)*deigv(2,:);
                            av=[av;av2]; dav=[dav;dav2];
                        else
                            av2=binom(eigv, -flg1+1); av2=av2'; av=[av;av2];
                        end
                end
        end
%
C=C_EBOF(flg1,flg2); D=[0 1 1;1 0 1;1 1 0];
%
if any(flg==[-3.2,-3.3])
    ut=1:10; vt=11:16; 
    u_ex=C(:,ut)*av(ut); v_ex=C(:,vt)*av(vt); A_ex=u_ex./v_ex;
else
    A_ex=C*av;
end
A_dv=D\(eigv-A_ex); A=diag([A_ex;A_dv],0);
for m=1:3, for n=1:3, if m~=n, k=9-m-n; A(m,n)=A(k,k); end , end, end
%
if nargout>1
    if any(flg==[-3.2,-3.3])
        du_ex=C(:,ut)*dav(ut,:); dv_ex=C(:,vt)*dav(vt,:);
        dA_ex=(v_ex.*du_ex-u_ex.*dv_ex)./v_ex.^2;
    else
        dA_ex=C*dav;
    end
    dA_dv=D\(deigv-dA_ex); dA_d=[dA_ex;dA_dv];
    for k=1:s3
        for i=1:6
            for j=1:6
                if i==j
                    dA6(i,j,k)=dA_d(i,k);
                else
                    if i<=3 && j<=3
                        dA6(i,j,k)=dA_d(9-i-j,k);
                    end
                end
            end
        end
    end
end
%
fnc=@(m,n) (m==n)*m+(m~=n)*(9-m-n);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                p1=fnc(i,j); p2=fnc(k,l);
                F(i,j,k,l)=A(p1,p2);
                if nargout>1, dF(i,j,k,l,:)=dA6(p1,p2,:); end
            end
        end
    end
end
%
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                ijkl=[i j k l]; val=0; 
                if nargout>1, dval=zeros(1,s3); end
                for r=1:3
                    for s=1:3
                        for m=1:3
                            for n=1:3
                                rsmn=[r s m n]; Q4=1; 
                                if nargout>1, dQ4=zeros(1,1,s3); end
                                for p=1:4
                                    Q4=Q4*Q(ijkl(p),rsmn(p)); 
                                    if nargout>1
                                        dQ4p=dQ(ijkl(p),rsmn(p),:);
                                        for q=1:4
                                            if p~=q
                                                dQ4p=dQ4p*Q(ijkl(q),rsmn(q));
                                            end
                                        end
                                        dQ4=dQ4+dQ4p;
                                    end
                                end
                                val=val+Q4*F(r,s,m,n);
                                if nargout>1
                                    for h=1:s3
                                        dval(h)=dval(h)+Q4*dF(r,s,m,n,h)+dQ4(h)*F(r,s,m,n);
                                    end
                                end
                            end
                        end
                    end
                end
                F2(i,j,k,l)=val;  if nargout>1, dF2(i,j,k,l,:)=dval;end
            end
        end
    end
end
end
end
%%
function [A4, dA4]=A4_IBOF(A,I,dA,dI,varargin)
A2v=v2M(A,2); d=I; s3=size(dA,3);
%
nA4=[1 1 1 1;1 1 2 2;1 1 2 3 ;...
     1 1 1 3;1 1 1 2;2 2 2 2 ;...
     2 2 2 3;2 2 1 3;2 2 1 2];

n4=size(nA4); A4v=zeros(n4(1),1); nn1=factorial(n4(2));
if nargout>1
    [bta, dbta]=bta_fnc(A,I,dA,dI,varargin{:}); dA4v=zeros(n4(1),s3);
else
    bta=bta_fnc(A,I,[],[],varargin{:});
end

for i=1:n4(1)
    ijkl=nA4(i,:);    A4v(i)=0;
    for j=1:3
        np=[perms(1:n4(2)) n4(2)+(1:j-1).*ones(nn1,1)];
        %
        nn=size(np); nm=nn(2)-j+1;
        np2=zeros(nn(1),2*(j+1)); np2(:,[1 end-nm+2:end])=np(:,1:nm);
        for jj=1:j-1
            jq1=2*j-4*jj+4;  jq2=nm-2*jj; j3=1/2*(j-1)*(j-2)*(jj-1);
            np2(:,jq1  :jq1+1)=repmat(np(:,nm+j-jj),1,2);
            np2(:,jq1+j3-2:jq1-1)=np(:,jq2+j3:jq2+1);
        end
        %
        ijklx=catj(ijkl ,j-1); nq=size(ijklx,1);
        %
        f0=0; df0=0; M2=repmat(A,1,1,j+1);
        for k=1:4-j
            f1=0; df1=0; kj=4*j-j*(j+1)/2+1-k;
            for m=1:nn(1)
                mn4=ijklx(:,np2(m,:)); 
                for r=1:nq
                    f2=1;
                    for n=1:j+1
                        p=2*n-1; 
                        f2=f2*M2(mn4(r,p),mn4(r,p+1),n);
                    end
                    f1=f1+f2;
                    %
                    if nargout>1
                        df2=0;
                        for n=1:j-k+2
                            p=2*n-1;
                            df3(1,:)=dA(mn4(r,2*j-p+2),mn4(r,2*j-p+3),:);
                            for s=1:j+1
                                if s~=n
                                    q=2*s-1;
                                    df3=df3*M2(mn4(r,2*j-q+2),mn4(r,2*j-q+3),j+2-s);
                                end
                            end
                            df2=df2+df3;
                        end
                        df1=df1+df2;
                    end
                end
            end
            f1= f1/nn(1) ;  f0=f0+bta(kj)*f1;  M2(:,:,k)=d;
            if nargout>1
                df1=df1/nn(1); df0=df0+bta(kj)*df1+dbta(kj,:)*f1;
            end
        end
        A4v(i)=A4v(i)+f0; if nargout>1, dA4v(i,:)=dA4v(i,:)+df0; end
    end
end
if nargout>1
    [A4,dA4]=A4v2M(A2v,A4v,dA4v,dA);
else
    A4=A4v2M(A2v,A4v);
end
end
%%
function [b,db]=bta_fnc(A,M,dA,dM,varargin)
s3=size(dA,3);
if nargout>1
    [I,dI]=INV(A,M,dA,dM,varargin{:});
else
    I=INV(A,M,[],[],varargin{:});
end
% I=[trace(A2);.5*(trace(A2)^2-trace(A2^2));det(A2)];
pv=[1,2,4,3,6,5,9,7,8,10,14,11,12,13,15,20,16,17,18,19,21];
b0=zeros(3,1); db0=zeros(3,s3);
for k=1:3
    b0(k)=C_IBOF(1,k);
    for i=1:5
        for j=0:i
            ij=j+.5*(i+1)*i+1; fac=C_IBOF(pv(ij),k);
            b0(k)=b0(k)+fac*I(2)^(i-j)*I(3)^j;
            if nargout>1
            db0(k,:)=db0(k,:)+fac*(i-j)*I(2)^(i-j-1)*I(3)^j*dI(2,:)+...
                              fac*j*I(2)^(i-j)*I(3)^(j-1)*dI(3,:);
            end
        end
    end
end

%
b=zeros(6,1); b([3 4 6])=b0;
b(1)=3/5*(-1/7+1/5*b(3)*(1/7+4/7*I(2)+8/3*I(3))-b(4)*(1/5-8/15*I(2)-14/15*I(3))+...
              -b(6)*(1/35-24/105*I(3)-4/35*I(2)+16/15*I(2)*I(3)+8/35*I(2)^2));
b(2)=6/7*(1-1/5*b(3)*(1+4*I(2))+7/5*b(4)*(1/6-I(2))-b(6)*(-1/5+2/3*I(3)+4/5*I(2)-8/5*I(2)^2));
b(5)=-4/5*b(3)-7/5*b(4)-6/5*b(6)*(1-4/3*I(2));
%
if nargout>1
db=zeros(6,s3); db([3 4 6],:)=db0;
db(1,:)=3/5*(1/5*db(3,:)*(1/7+4/7*I(2)+8/3*I(3))-db(4,:)*(1/5-8/15*I(2)-14/15*I(3))+...
                -db(6,:)*(1/35-24/105*I(3)-4/35*I(2)+16/15*I(2)*I(3)+8/35*I(2)^2)+...
             1/5*b(3)*(4/7*dI(2,:)+8/3*dI(3,:))-b(4)*(-8/15*dI(2,:)-14/15*dI(3,:))+...
                 -b(6)*(-24/105*dI(3,:)-4/35*dI(2,:)+16/15*(dI(2,:)*I(3)+I(2)*dI(3,:))+...
                 16/35*I(2)*dI(2,:)));
db(2,:)=6/7*(-1/5*db(3,:)*(1+4*I(2))+7/5*db(4,:)*(1/6-I(2))-db(6,:)*(-1/5+2/3*I(3)+4/5*I(2)-8/5*I(2)^2)+...
             -4/5*b(3)*dI(2,:)-7/5*b(4)*dI(2,:)-b(6)*(2/3*dI(3,:)+4/5*dI(2,:)-16/5*I(2)*dI(2,:)));
db(5,:)=-4/5*db(3,:)-7/5*db(4,:)-6/5*db(6,:)*(1-4/3*I(2))+8/5*b(6)*dI(2,:);
end
end
%%
function [f,df]=INV(A,I,dA,dI,varargin)
if nargout>1
    [v, ~, dv]=deigvdA(A,I,dA,dI,varargin{:});
else
    v=deigvdA(A,I,[],[],varargin{:});
end
f=zeros(3,1);
for i=1:3
    ij=nchoosek(1:3,i); nij=size(ij);
    f0=0; df0=0;
    for j=1:nij(1)
        f1=1; df1=0;
        for k=1:nij(2)
            f1=f1*v(ij(j,k));
            %
            if nargout>1
                df2=dv(ij(j,k),:);
                for m=1:nij(2)
                    if m~=k
                        df2=df2*v(ij(j,m));
                    end
                end
                df1=df1+df2;
            end
            %
        end
        f0=f0+f1; if nargout>1,  df0=df0+df1; end
    end
    f(i)=f0; if nargout>1, df(i,:)=df0; end
end
end
%%
function [f,df]=A4v2M(A2v,A4v,varargin)
% A2v=['A11','A12','A13','A22','A23']
% A4v=['A1111','A1122','A1123','A1113','A1112','A2222','A2223','A2213','A2212']
d=eye(3); f1=@(i,j) (9-i-j).*(1-d(i+3*(j-1)))+d(i+3*(j-1)).*i;
f2=@(m) [m, m]+(m>3)*([-m, 9-2*m]+[1 -1].*(m^2-11*m+32)/2);
Mij=@(M,i,j) logical(prod((M==i)+(M==j),2)); 
%
jj=0; Axx=[];B2=[];
for i =1:3
    for j=i:3
        for k=j:3
            for l=k:3
                ijkl=[i j k l]; jj=jj+1;  Axx=[Axx; f1(i,j) f1(k,l)];
                B1=unique(perms(ijkl),'Rows'); njj=jj*ones(size(B1,1),1);
                B2=[B2;njj sum((3.^(0:3)).*(B1-1),2)+1];
            end
        end
    end
end
%
k=0; for i=1:3,for j=i:6, k=k+1; D1(i,j)=k; end, end
%
Axx(Mij(Axx,4,6),:)=[2 5]; Axx=sort(Axx,2);
[~,ndk]=sortrows(Axx,[1 2]); D2=(1:15)'; D3=D2(ndk,1);
%
jj=0;
for i=1:2
    for j=i:6
        if j~=3
            jj=jj+1;
            n_inp(jj,1)=D3(D1(i, j),1);
        end
    end
end
D2(n_inp,2)=A4v; A2=v2M(A2v,1);
if nargin>2
    [dA4v,dA]=varargin{1:2}; dD2(n_inp,:)=dA4v; s3=size(dA,3); 
end
%
for k=1:2
    for i=1:3
        m=i+3*(k-1); r=f2(m); s=[(2-k) (k-1)]*[m 3;3 m];
        val=A2(r(1),r(2));
        if nargin>2, dval(1,:,1)=dA(r(1),r(2),:); end
        for j=1:2
            p=sort([j m],2);
            val=val-D2(D3(D1(p(1),p(2)),1),2);
            if nargin>2
                dval=dval-dD2(D3(D1(p(1),p(2)),1),:);
            end
        end
        D2(D3(D1(s(1),s(2)),1),2)=val;
        if nargin>2
            dD2(D3(D1(s(1),s(2)),1),:)=dval;
        end
    end
end
f(B2(:,2))= D2(B2(:,1),2); f=reshape(f,3,3,3,3);
%
if nargout>1
    df(B2(:,2),:)= dD2(B2(:,1),:); df=reshape(df,3,3,3,3,s3);
end
%
end
%%
function [v,Q, dv,dQ]=deigvdA(K,M,dK,dM,k,flg,varargin)
switch flg
    case {1 2 3}
        % syms v; v=vpasolve(det(K-v*M)==0,v); v=double(sort(v,'descend'));
        % det(K)-sum(M.*adjoint(K)','all')*v+sum(K.*adjoint(M)','all')*v^2-det(M)*v^3=0;
        v=roots([-det(M),sum(K.*adjoint(M)','all'),-sum(M.*adjoint(K)','all'),det(K)]);
        v=sort(v,'descend');
        %
        [s1,s2]=size(K); if nargout>2,s3=size(dK,3); end
        for i=1:s2
            S=K-v(i)*M;
            Qi=rand(s1,1); err=1;
            while err>1e-8
                Ci=Qi'*M*Qi;
                switch flg
                    case 1
                        Rk=Ci-1; dSk=2*Qi'*M;
                    case 2
                        [p,a]=varargin{1:2};Z=zeros(1,s2); Z(p)=1;
                        Rk=Z*Qi-a; dSk=Z;
                    case 3
                        b=varargin{1};
                        Rk=sqrt(Qi'*Qi)-b; dSk=Qi';
                end
                %
                R=S*Qi; R(k)=Rk; J=S; J(k,:)=dSk;
                Qn=Qi-2.*J\R; err=norm(Qn-Qi); Qi=Qn;
            end
            Q(:,i)=Qi;
            %
            if nargout>2
                for j=1:s3
                    dv(i,j)=1/Ci*(Qi'*(dK(:,:,j)-v(i)*dM(:,:,j))*Qi);
                    F=-(dK(:,:,j)-dv(i,j)*M-v(i)*dM(:,:,j))*Qi;
                    switch flg
                        case 1
                            Fk=-Qi'*dM(:,:,j)*Qi;
                        case {2 3}
                            Fk=0;
                    end
                    dF=F; dF(k)=Fk; P=S; P(k,:)=dSk; dQ(:,i,j)=P\dF;
                end
            end
        end
    case 4
        if nargout>2
            [v, Q, dv, dQ]=Nlsn(K,M,dK,dM,k);
        else
            [v, Q]=Nlsn(K,M,dK,dM,k);
        end
end
end
%%
function [v, Q, dv, dQ]=Nlsn(K,M,dK,dM,k)
[Q, v]=eig(K);   [v, ind]=sort(diag(v,0),'descend');Q=Q(:,ind);
[s1,s2]=size(K); 
if nargout>2
    s3=size(dK,3); dv=zeros(s2,s3); dQ=zeros(s1,s2,s3);
    for j=1:s2
        L=v(j); Qj=Q(:,j);
        Cj=Qj'*M*Qj; Sj=K-L*M; Sj(k,:)=0; Sj(:,k)=0; Sj(k,k)=1;
        for i=1:s3
            dKi=dK(:,:,i); dMi=dM(:,:,i);
            dLi=1/Cj*(Qj'*(dKi-L*dMi)*Qj); dv(j,i)=dLi;
            if nargout>3
                Fj=-(dKi-dLi*M-L*dMi)*Qj; Fj(k)=0;
                Vj=Sj\Fj; cj=-Qj'*M*Vj-1/2*Qj'*dMi*Qj;
                dQ(:,j,i)=Vj+cj.*Qj;
            end
        end
    end
end
end
%%
function [detj,ddetj]=f_detj(A,dA)
detj=0; s1=size(A,1); 
if nargout> 1, s3=size(dA,3); ddetj=zeros(1,1,s3); end
ijkn=perms(1:s1); [s1, s2]=size(ijkn);
for i=1:s1
    f=1; ijk=ijkn(i,:); eijk=par(ijk);
    for m=1:s2
        f=f*A(ijk(m),m);
        if nargout>1
            df=dA(ijk(m),m,:);
            for n=1:s2
                if m~=n
                    df=df*A(ijk(n),n);
                end
            end
            ddetj=ddetj+eijk*df;
        end
    end
    detj=detj+eijk*f;
end
end
%%
function [R,J]=dMdA(func, A, del, flg, varargin)
A=A(:); n2=numel(A); nA=size(A); nA=num2cell(nA);
D=zeros(n2,1); tol=eps;
f=@(A) func(A,varargin{:}); R=f(A); n1=numel(R);
%
J=zeros(n1*n2,1);
Fac={ 
    1,         [-1  1 ; -1  0] ;
    1,         [-1  1 ;  0  1] ;
    2,         [-1  1 ; -1  1] ;
    2,     [1 -4  3 ; -2 -1 0] ;
    2,     [-3 4 -1 ;  0  1 2] ;
   12, [1 -8  8 -1 ; -2 -1 1 2]
                                    };
%
M1=Fac{flg,1}; M2=Fac{flg,2}; m=length(M2); 
for i=1:n2
    na=(i-1)*n1+1; nb=na+n1-1;
    D(i,1)=1; J(na:nb,1)=zeros(n1,1); dA=del*D; h=del; 
    if abs(A(i))>tol,  dA=dA.*A(:); h=h*A(i);  end
    for j=1:m
        Jj=M2(1,j)*f(A+M2(2,j)*reshape(dA,nA{:})); 
        J(na:nb,1)=J(na:nb,1)+Jj(:);
    end
    J(na:nb,1)=J(na:nb,1)/(M1*h);
    D(i,1)=0;
end
nR=size(R); nR=num2cell(nR);J=reshape(J,nR{:},n2);
end
%%
function [f,df]=gencls(A,dA,flg)
d=eye(3); f=zeros(3,3,3,3);
%
if nargout>1
    s3=size(dA,3); [b,db]=fbta(A,dA,flg); df=zeros(3,3,3,3,s3);
else
    b=fbta(A,[],flg); 
end
%
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                f(i,j,k,l)=b(1)*(d(i,j)*d(k,l))+b(2)*(d(i,k)*d(j,l)+d(i,l)*d(j,k))+...
                    b(3)*(d(i,j)*A(k,l)+A(i,j)*d(k,l))+...
                    b(4)*(A(i,k)*d(j,l)+A(j,l)*d(i,k)+A(i,l)*d(j,k)+A(j,k)*d(i,l))+...
                    b(5)*(A(i,j)*A(k,l))+b(6)*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
                if nargout>1
                    df(i,j,k,l,:)=db(1,1,:)*(d(i,j)*d(k,l))+db(1,2,:)*(d(i,k)*d(j,l)+d(i,l)*d(j,k))+...
                        db(1,3,:)*(d(i,j)*A(k,l)+A(i,j)*d(k,l))+...
                        db(1,4,:)*(A(i,k)*d(j,l)+A(j,l)*d(i,k)+A(i,l)*d(j,k)+A(j,k)*d(i,l))+...
                        db(1,5,:)*(A(i,j)*A(k,l))+db(1,6,:)*(A(i,k)*A(j,l)+A(i,l)*A(j,k));
                    for q=1:s3
                        df(i,j,k,l,q)=df(i,j,k,l,q)+...
                            b(3)*(d(i,j)*dA(k,l,q)+dA(i,j,q)*d(k,l))+...
                            b(4)*(dA(i,k,q)*d(j,l)+dA(j,l,q)*d(i,k)+dA(i,l,q)*d(j,k)+dA(j,k,q)*d(i,l))+...
                            b(5)*(dA(i,j,q)*A(k,l)+A(i,j)*dA(k,l,q))+...
                            b(6)*(dA(i,k,q)*A(j,l)+A(i,k)*dA(j,l,q)+dA(i,l,q)*A(j,k)+A(i,l)*dA(j,k,q));
                    end
                end
                %
                for m=1:3
                    f(i,j,k,l)=f(i,j,k,l)+b(7)*(d(i,j)*A(k,m)*A(m,l)+A(i,m)*A(m,j)*d(k,l));
                    if nargout>1
                        for q=1:s3
                            df(i,j,k,l,q)=df(i,j,k,l,q)+db(1,7,q)*(d(i,j)*A(k,m)*A(m,l)+A(i,m)*A(m,j)*d(k,l))+...
                                b(7)*(d(i,j)*dA(k,m,q)*A(m,l)+d(i,j)*A(k,m)*dA(m,l,q)+...
                                dA(i,m,q)*A(m,j)*d(k,l)+A(i,m)*dA(m,j,q)*d(k,l));
                        end
                    end
                    for n=1:3
                        f(i,j,k,l)=f(i,j,k,l)+b(8)*(A(i,m)*A(m,j)*A(k,n)*A(n,l));
                        if nargout>1
                            for q=1:s3
                                df(i,j,k,l,q)=df(i,j,k,l,q)+db(1,8,q)*(A(i,m)*A(m,j)*A(k,n)*A(n,l))+...
                                    b(8)*(dA(i,m,q)*A(m,j)*A(k,n)*A(n,l)+A(i,m)*dA(m,j,q)*A(k,n)*A(n,l)+...
                                    A(i,m)*A(m,j)*dA(k,n,q)*A(n,l)+A(i,m)*A(m,j)*A(k,n)*dA(n,l,q));
                            end
                        end
                    end
                end
                %
            end
        end
    end
end
end
%%
function [b,db]=fbta(A,dA,flg)
%
ysq=sum(A.*A','all'); a=exp(2*(1-3*ysq)/(1-ysq));
%
bta=[
    1/15        1/15         0           0       0      0        0        0;
   -1/35       -1/35       1/7         1/7       0      0        0        0;
       0           0         0           0       1      0        0        0;
       0           0         0           0       1      1        0   -2/ysq;
       0           0       2/5           0    -1/5    3/5     -2/5        0;
26*a/315    26*a/315   16*a/63     -4*a/21       1      1        0   -2/ysq];
%
b=bta(flg,:);
%
if nargout>1
    s3=size(dA,3); dysq=sum(2*dA.*A',[1 2]); da=(-4*a/(1-ysq)^2)*dysq; 
    dbta=zeros(6,8,s3);  dbta([4 6],8,:)=2/ysq^2*dysq.*[1;1];
    dbta(6,1:4,:)=[26/315    26/315   16/63    -4/21].*da; db=dbta(flg,:,:);
end
end
%%
function f=DA(d)
for i=1:3
    for j=1:3
        for m=1:2
            for n=m:3
                k=2*(m-1)+n; p=-.5*(i^2-3*i+2);
                f(i,j,k)=d(i,m)*d(j,n)+(1-d(i,j))*d(i,n)*d(j,m)+p*d(i,j)*d(m,n);
            end
        end
    end
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
%%
function A=catj(A,p)
pj=1;
while pj<=p
    n=size(A); B=zeros(3*n(1),n(2)+1);
    for k=1:3
        rk=(k-1)*n(1)+1:k*n(1);
        for j=1:n(1)
            B(rk,1:n(2))=A; B(rk,n(2)+1)=k;
        end
    end
    A=B; pj=pj+1;
end
end
%%
function [t,X]=rk4(fnc,t, X0, varargin)
n=length(t);  X(1,:)=X0';
for j=1:n-1
    dt=t(j+1)-t(j); 
    K1=fnc(t(j)     , X(j,:)'         ,varargin{:});
    K2=fnc(t(j)+dt/2, X(j,:)' +dt/2*K1,varargin{:});
    K3=fnc(t(j)+dt/2, X(j,:)' +dt/2*K2,varargin{:});
    K4=fnc(t(j)+dt  , X(j,:)' +dt*K3  ,varargin{:});
    dXj=dt/6*(K1+2*K2+2*K3+K4); X(j+1,:) = X(j,:) + dXj';
end
end
%%
function [av,dav]=binom(v, n)
av=1; if nargout>1, dav=[0;0]; end
for i=1:n
    for j=0:i
        av=[av v(1)^(i-j)*v(2)^j];
        if nargout>1
            dav=[dav,...
                [(i-j)*v(1)^(i-j-1)*v(2)^j      ;
                j*v(1)^(i-j  )*v(2)^(j-1)]];
        end
    end
end
end
%%
function [f,a]=par(a)
n=length(a); k=0; j=1; m=0;
while 1
    aj=a(j); 
    if j~=aj
        c=a(aj); a(aj)=aj;  a(j)=c; k=k+1; m=m+abs(a(j)-j);
    end
    m=m*(j<n); j=j*(j<n)+1; 
    if (j==n)*(m==0)==1, break, end
end
f=(-1)^k;
end
%%
function f=C_IBOF(i,j)
C=[ 
 0.249409081657860E+02  -0.497217790110754E+00   0.234146291570999E+02
-0.435101153160329E+03   0.234980797511405E+02  -0.412048043372534E+03
 0.372389335663877E+04  -0.391044251397838E+03   0.319553200392089E+04
 0.703443657916476E+04   0.153965820593506E+03   0.573259594331015E+04
 0.823995187366106E+06   0.152772950743819E+06  -0.485212803064813E+05
-0.133931929894245E+06  -0.213755248785646E+04  -0.605006113515592E+05
 0.880683515327916E+06  -0.400138947092812E+04  -0.477173740017567E+05
-0.991630690741981E+07  -0.185949305922308E+07   0.599066486689836E+07
-0.159392396237307E+05   0.296004865275814E+04  -0.110656935176569E+05
 0.800970026849796E+07   0.247717810054366E+07  -0.460543580680696E+08
-0.237010458689252E+07   0.101013983339062E+06   0.203042960322874E+07
 0.379010599355267E+08   0.732341494213578E+07  -0.556606156734835E+08
-0.337010820273821E+08  -0.147919027644202E+08   0.567424911007837E+09
 0.322219416256417E+05  -0.104092072189767E+05   0.128967058686204E+05
-0.257258805870567E+09  -0.635149929624336E+08  -0.152752854956514E+10
 0.214419090344474E+07  -0.247435106210237E+06  -0.499321746092534E+07
-0.449275591851490E+08  -0.902980378929272E+07   0.132124828143333E+09
-0.213133920223355E+08   0.724969796807399E+07  -0.162359994620983E+10
 0.157076702372204E+10   0.487093452892595E+09   0.792526849882218E+10
-0.232153488525298E+05   0.138088690964946E+05   0.466767581292985E+04
-0.395769398304473E+10  -0.160162178614234E+10  -0.128050778279459E+11];
f=C(i,j);
end
%%
function C=C_EBOF(flg1,flg2)
switch flg1
    case -1 % linear
        switch flg2
            case 1 % Smooth ortho - Cintra and Tucker
                C =[-0.15  1.15 -0.10 ;
                    -0.15  0.15  0.90 ;
                     0.60 -0.60 -0.60];
            case 2 % General ortho - Kuzmin
                C=1/7*[ -3/5  6   0;
                        -3/5  0   6;
                        27/5 -6  -6];
        end
    case -2 % quadratic 
        %
        switch flg2
            case 1 % fitted ortho - Cintra & Tucker
                C=[
                    0.060964  0.371243 -0.369160 0.555301 0.371218 0.318266 ;
                    0.124711 -0.389402  0.086169 0.258844 0.544992 0.796080 ;
                    1.228982 -2.054116 -2.260574 0.821548 1.819756 1.053907];
            case 2 % general ortho - Cintra and Tucker
                C=[0  0  0 1  0  0 ;
                   0  0  0 0  0  1 ;
                   1 -2 -2 1  2  1];
            case 3 % natural ortho - exact midpoint fit -Kuzmin
                C=[0.0708  0.3236  -0.3776 0.6056  0.4124  0.3068;
                   0.0708 -0.2792   0.2252 0.2084  0.4124  0.7040;
                   1.1880 -2.0136  -2.1264 0.8256  1.7640  0.9384];
            case 4 % Improved ortho fitted ORW Chung et al 2002
                C=[0.070055  0.339376 -0.396796 0.590331 0.411944 0.333693 ;
                   0.115177 -0.368267  0.094820 0.252880 0.535224 0.800181 ;
                   1.249811 -2.148297 -2.290157 0.898521 1.934914 1.044147];
        end
    case -3 % Cubic
        switch flg2
            case 1 % Natural ortho - extended quad. fit
                C=[0  1/2     0 1/2  -3/5   0  0  3/5 3/5 0;
                   0    0   1/2   0  -3/5 1/2  0  3/5 3/5 0;
                   1 -3/2  -3/2 1/2   2/5 1/2  0  3/5 3/5 0];
            case {2 3} % Rational ellipsoid 
                switch flg2
                    case 2 % Wertzel
                        C=[
                             0.1433751825  0.1433751825  0.9685744898
                            -0.6566650339 -0.5209453949 -2.5526857671
                            -0.5106016916 -0.6463213306 -2.5756669706
                             4.4349137241  2.3303190917  4.4520903005
                             3.5295952199  0.6031924921  2.2044050704
                             0.1229618909  5.1539592511  2.2485545147
                            -5.5556896198 -1.6481269200 -1.8811803355
                            -2.8284365891 -5.4494528976 -1.9023485762
                            -2.9144388828 -0.2256222796 -0.6202937932
                             0.2292109036 -3.7461520908 -0.6414620339
                             1.0000000000  1.0000000000  1.0000000000
                             0.7257989503  0.6916858207 -1.2134964928
                             3.0941511876  3.1282643172 -1.2128608265
                            -4.7303686308 -4.7303686308  0.6004510415
                            -1.6239324646 -1.5898193351  0.2393747647
                            -3.1742364608 -3.2083495904  0.2162486576]';        
                    case 3 %LAR32 Matthew Mullens
                        C=[ 
                             0.0876022	 0.1568052	 1.0724237
                             0.0282056	-0.5778189	-2.8035540
                            -0.4267843	-0.5142809	-2.6615761
                             0.8764691	 2.1323050	 4.5667285
                             1.2746771	 0.6842509	 2.3893798
                             0.6020316	 3.4548353	 2.0975231
                            -1.9189311	-1.6141226	-1.9047047
                            -0.9342913	-4.0052611	-1.7549784
                            -1.0665831	-0.2632371	-0.6582489
                            -0.2628549	-2.2281332	-0.5082827
                             1.0000000	 1.0000000	 1.0000000
                            -0.2440019	 0.3656529	-1.0685125
                            -0.5741509	 1.3857255	-0.7713565
                            -0.8952261	-2.8663578	 0.2069083
                            -0.4320974	-1.3596872	 0.0673869
                            -0.4627095	-1.5189962	-0.2489999]';

                end
                pt=[1 2 3 5 4 6 9 7 8 10 11 12 13 15 14 16];
                C=C(:,pt);
            case 4 % Improved Orthotropic fitted - ORW3 - Chung & kwon
                C=[ -0.1480648093 -0.2106349673  0.4868019601
                     0.8084618453  0.9092350296  0.5776328438
                     0.3722003446 -1.2840654776 -2.2462007509
                     0.7765597096  1.1104441966  0.4605743789
                    -1.3431772379  0.1260059291	-1.9088154281
                    -1.7366749542 -2.5375632310	-4.8900459209
                     0.8895946393  1.9988098293	 4.0544348937
                     1.7367571741  1.4863151577	 3.8542602127
                    -0.0324756095  0.5856304774	 1.1817992322
                     0.6631716575 -0.0756740034	 0.9512305286]';
                 pt=[1 2 4 3 6 5 9 7 8 10]; C=C(:,pt);
        end
    case -4 % quartic
        switch flg2
            case 1 % Verweyst
                C=[
                    0.6363    0.6363    2.7405;
                   -1.8727   -3.3153   -9.1220;
                   -4.4797   -3.0371  -12.2571;
                   11.9590   11.8273   34.3199;
                    3.8446    6.8815   13.8295;
                   11.3421    8.4368   25.8685;
                  -10.9583  -15.9121  -37.7029;
                  -20.7278  -15.1516  -50.2756;
                   -2.1162   -6.4873  -10.8802;
                  -12.3876   -8.6389  -26.9637;
                    9.8160    9.3252   27.3347;
                    3.4790    7.7468   15.2651;
                   11.7493    7.4815   26.1135;
                    0.5080    2.2848    3.4321;
                    4.8837    3.5977   10.6117]';
            case 2 % FFLAR4 - Matthew Mullens
                C=[ 
                      0.678225884   0.748226727   3.167356369
                     -3.834359034  -4.249612053 -13.288266400
                     -2.664862865  -2.987266447 -11.680179330
                     14.209962670  14.938209410  43.700607680
                      9.746185193   8.641488072  23.788431340
                      2.700369681   5.974489008  17.383121430
                    -22.447252700 -21.757217160 -58.354308000
                    -13.078649640 -15.798676320 -49.513705640
                     -8.013024236  -7.521216405 -19.959054610
                     -0.125467689  -3.616551654 -11.755525930
                     12.689484570  12.640352670  35.425354130
                     10.563248410  10.222185780  25.844317920
                      2.487386515   4.788201652  18.226443930
                      2.417857515   2.376441613   6.291273472
                     -0.328195677   1.056519961   2.925785795]';
            case 3  %LAR4 Ref Matthew Mullens
                C=[  
                      0.813175172   1.768619587   4.525066937
                     -3.065410883  -9.826017151 -19.259137620
                     -4.659333003  -6.484058476 -17.650178090
                     14.747639770  28.905936750  61.543979540
                      6.329870878  19.986994700  33.901239610
                      9.739797775  10.759963010  28.467355970
                    -15.922240910 -40.492387100 -76.738638810
                    -20.818571900 -27.442217500 -68.977583290
                     -4.216519964 -17.715409270 -27.768082700
                     -8.993993112  -7.230748101 -22.399036130
                     11.470974520  19.729631240  43.875135630
                      5.834142985  18.709047480  32.480679940
                      9.874209286   8.882877701  26.928320210
                      1.138888034   5.785725498	  8.600822308
                      3.100457733   2.224834058	  7.101978254]';

        end
        pt=[1 2 3 5 4 6 9 7 8 10 14 12 11 13 15];
        C=C(:,pt);
end
end
%