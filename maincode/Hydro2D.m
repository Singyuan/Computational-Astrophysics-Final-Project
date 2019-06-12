clc; clear all;

% Strong Shock Tube Problem, Flash 7.2
rhoL=1; uL=0; vL=0; PL=1; rhoR=0.125; uR=0; vR=0; PR=0.1; shift=50; M=100; dx=1/(100-1); dy=1/(100-1); gamma=1.4;
tboundary=0.2; 

rho=zeros(M,M); u=zeros(M,M); v=zeros(M,M); P=zeros(M,M); e=zeros(M,M);
for i=1:M
for j=1:M
    if j>i
        rho(j,i)=rhoL;
        u(j,i)=uL;
        v(j,i)=vL;
        P(j,i)=PL;
    else
        rho(j,i)=rhoR;
        u(j,i)=uR;
        v(j,i)=vR;
        P(j,i)=PR;
    end
end
end
for i=1:M
for j=1:M
    e(j,i)=P(j,i)/(gamma-1)/rho(j,i);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhomarch=zeros(M,M); umarch=zeros(M,M); vmarch=zeros(M,M); Pmarch=zeros(M,M); emarch=zeros(M,M);
rho1=zeros(M,M); u1=zeros(M,M); v1=zeros(M,M); P1=zeros(M,M); e1=zeros(M,M);
t=0;
while t<tboundary
    Smax=0;
    for i=2:M
        for j=2:M
            [SL1,SR1,Sstar1,FL1,FR1,FstarL1,FstarR1] = VectorsX(rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),v(i,j-1),v(i,j),P(i,j-1),P(i,j));
            [SL3,SR3,Sstar3,GL1,GR1,GstarL1,GstarR1] = VectorsY(rho(i-1,j),rho(i,j),u(i-1,j),u(i,j),v(i-1,j),v(i,j),P(i-1,j),P(i,j));
            if Smax<max([abs(SL1) abs(SR1) abs(SL3) abs(SR3)])
                Smax=max([abs(SL1) abs(SR1) abs(SL3) abs(SR3)]);
            end
        end
    end
    Ccfl=0.5; 
    d=min(dx,dy); dt=Ccfl*d/Smax;
    
    for i=1:M
        for j=2:(M-1)  
            [SL1,SR1,Sstar1,FL1,FR1,FstarL1,FstarR1] = VectorsX(rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),v(i,j-1),v(i,j),P(i,j-1),P(i,j));
            [SL2,SR2,Sstar2,FL2,FR2,FstarL2,FstarR2] = VectorsX(rho(i,j),rho(i,j+1),u(i,j),u(i,j+1),v(i,j),v(i,j+1),P(i,j),P(i,j+1));

            E=rho(i,j)*(u(i,j).^2/2+v(i,j).^2/2+e(i,j));
            U=[rho(i,j) rho(i,j)*u(i,j) rho(i,j)*v(i,j) E];
            Fminus=HLLCFlux(SL1,SR1,Sstar1,FL1,FR1,FstarL1,FstarR1);
            Fplus=HLLCFlux(SL2,SR2,Sstar2,FL2,FR2,FstarL2,FstarR2);
            U=U+dt/dx*(Fminus-Fplus);
            rhomarch(i,j)=U(1,1);
            umarch(i,j)=U(1,2)./rhomarch(i,j);
            vmarch(i,j)=U(1,3)./rhomarch(i,j);
            emarch(i,j)=U(1,4)./rhomarch(i,j)-umarch(i,j).^2/2-vmarch(i,j).^2/2;
            Pmarch(i,j)=emarch(i,j)*(gamma-1).*rhomarch(i,j); 
        end
        rhomarch(i,1)=rhomarch(i,2); rhomarch(i,M)=rhomarch(i,M-1);
        umarch(i,1)=umarch(i,2); umarch(i,M)=umarch(i,M-1);
        vmarch(i,1)=vmarch(i,2); vmarch(i,M)=vmarch(i,M-1);
        emarch(i,1)=e(i,2); emarch(i,M)=emarch(i,M-1);
        Pmarch(i,1)=Pmarch(i,2); Pmarch(i,M)=Pmarch(i,M-1);
        rho1(i,:)=rhomarch(i,:);
        u1(i,:)=umarch(i,:);
        v1(i,:)=vmarch(i,:);
        P1(i,:)=Pmarch(i,:);
        e1(i,:)=emarch(i,:);
    end
    
    for j=1:M
        for i=2:(M-1)
            [SL3,SR3,Sstar3,GL1,GR1,GstarL1,GstarR1] = VectorsY(rho1(i-1,j),rho1(i,j),u1(i-1,j),u1(i,j),v1(i-1,j),v1(i,j),P1(i-1,j),P1(i,j));
            [SL4,SR4,Sstar4,GL2,GR2,GstarL2,GstarR2] = VectorsY(rho1(i,j),rho1(i+1,j),u1(i,j),u1(i+1,j),v1(i,j),v1(i+1,j),P1(i,j),P1(i+1,j));
            
            E=rho1(i,j)*(u1(i,j).^2/2+v1(i,j).^2/2+e1(i,j));
            U=[rho1(i,j) rho1(i,j)*u1(i,j) rho1(i,j)*v1(i,j) E];
            Gminus=HLLCFlux(SL3,SR3,Sstar3,GL1,GR1,GstarL1,GstarR1);
            Gplus=HLLCFlux(SL4,SR4,Sstar4,GL2,GR2,GstarL2,GstarR2);
            U=U+dt/dy*(Gminus-Gplus);
            rhomarch(i,j)=U(1,1);
            umarch(i,j)=U(1,2)./rhomarch(i,j);
            vmarch(i,j)=U(1,3)./rhomarch(i,j);
            emarch(i,j)=U(1,4)./rhomarch(i,j)-umarch(i,j).^2/2-vmarch(i,j).^2/2;
            Pmarch(i,j)=emarch(i,j)*(gamma-1).*rhomarch(i,j); 
        end
        rhomarch(1,j)=rhomarch(2,j); rhomarch(M,j)=rhomarch(M-1,j);
        umarch(1,j)=umarch(2,j); umarch(M,j)=umarch(M-1,j);
        vmarch(1,j)=vmarch(2,j); vmarch(M,j)=vmarch(M-1,j);
        emarch(1,j)=emarch(2,j); emarch(M,j)=emarch(M-1,j);
        Pmarch(1,j)=Pmarch(2,j); Pmarch(M,j)=Pmarch(M-1,j);
        rho1(:,j)=rhomarch(:,j);
        u1(:,j)=umarch(:,j);
        v1(:,j)=vmarch(:,j);
        P1(:,j)=Pmarch(:,j);
        e1(:,j)=emarch(:,j);
    end
    
    for i=1:M
        for j=2:(M-1)  
            [SL1,SR1,Sstar1,FL1,FR1,FstarL1,FstarR1] = VectorsX(rho1(i,j-1),rho1(i,j),u1(i,j-1),u1(i,j),v1(i,j-1),v1(i,j),P1(i,j-1),P1(i,j));
            [SL2,SR2,Sstar2,FL2,FR2,FstarL2,FstarR2] = VectorsX(rho1(i,j),rho1(i,j+1),u1(i,j),u1(i,j+1),v1(i,j),v1(i,j+1),P1(i,j),P1(i,j+1));
            [SL3,SR3,Sstar3,FL3,FR3,FstarL3,FstarR3] = VectorsX(rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),v(i,j-1),v(i,j),P(i,j-1),P(i,j));
            [SL4,SR4,Sstar4,FL4,FR4,FstarL4,FstarR4] = VectorsX(rho(i,j),rho(i,j+1),u(i,j),u(i,j+1),v(i,j),v(i,j+1),P(i,j),P(i,j+1));

            E=rho(i,j)*(u(i,j).^2/2+v(i,j).^2/2+e(i,j));
            U=[rho(i,j) rho(i,j)*u(i,j) rho(i,j)*v(i,j) E];
            Fminusprime=HLLCFlux(SL1,SR1,Sstar1,FL1,FR1,FstarL1,FstarR1);
            Fplusprime=HLLCFlux(SL2,SR2,Sstar2,FL2,FR2,FstarL2,FstarR2);
            Fminus=HLLCFlux(SL3,SR3,Sstar3,FL3,FR3,FstarL3,FstarR3);
            Fplus=HLLCFlux(SL4,SR4,Sstar4,FL4,FR4,FstarL4,FstarR4);
            U=U+dt/dx*((Fminus+Fminusprime)/2-(Fplus+Fplusprime)/2);
            rhomarch(i,j)=U(1,1);
            umarch(i,j)=U(1,2)./rhomarch(i,j);
            vmarch(i,j)=U(1,3)./rhomarch(i,j);
            emarch(i,j)=U(1,4)./rhomarch(i,j)-umarch(i,j).^2/2-vmarch(i,j).^2/2;
            Pmarch(i,j)=emarch(i,j)*(gamma-1).*rhomarch(i,j); 
        end
        rhomarch(i,1)=rhomarch(i,2); rhomarch(i,M)=rhomarch(i,M-1);
        umarch(i,1)=umarch(i,2); umarch(i,M)=umarch(i,M-1);
        vmarch(i,1)=vmarch(i,2); vmarch(i,M)=vmarch(i,M-1);
        emarch(i,1)=e(i,2); emarch(i,M)=emarch(i,M-1);
        Pmarch(i,1)=Pmarch(i,2); Pmarch(i,M)=Pmarch(i,M-1);
        rho(i,:)=rhomarch(i,:);
        u(i,:)=umarch(i,:);
        v(i,:)=vmarch(i,:);
        P(i,:)=Pmarch(i,:);
        e(i,:)=emarch(i,:);
    end
    
    for j=1:M
        for i=2:(M-1)
            [SL1,SR1,Sstar1,GL1,GR1,GstarL1,GstarR1] = VectorsY(rho1(i-1,j),rho1(i,j),u1(i-1,j),u1(i,j),v1(i-1,j),v1(i,j),P1(i-1,j),P1(i,j));
            [SL2,SR2,Sstar2,GL2,GR2,GstarL2,GstarR2] = VectorsY(rho1(i,j),rho1(i+1,j),u1(i,j),u1(i+1,j),v1(i,j),v1(i+1,j),P1(i,j),P1(i+1,j));
            [SL3,SR3,Sstar3,GL3,GR3,GstarL3,GstarR3] = VectorsY(rho(i-1,j),rho(i,j),u(i-1,j),u(i,j),v(i-1,j),v(i,j),P(i-1,j),P(i,j));
            [SL4,SR4,Sstar4,GL4,GR4,GstarL4,GstarR4] = VectorsY(rho(i,j),rho(i+1,j),u(i,j),u(i+1,j),v(i,j),v(i+1,j),P(i,j),P(i+1,j));
            
            E=rho(i,j)*(u(i,j).^2/2+v(i,j).^2/2+e(i,j));
            U=[rho(i,j) rho(i,j)*u(i,j) rho(i,j)*v(i,j) E];
            Gminusprime=HLLCFlux(SL1,SR1,Sstar1,GL1,GR1,GstarL1,GstarR1);
            Gplusprime=HLLCFlux(SL2,SR2,Sstar2,GL2,GR2,GstarL2,GstarR2);
            Gminus=HLLCFlux(SL3,SR3,Sstar3,GL3,GR3,GstarL3,GstarR3);
            Gplus=HLLCFlux(SL4,SR4,Sstar4,GL4,GR4,GstarL4,GstarR4);
            U=U+dt/dy*((Gminus+Gminusprime)/2-(Gplus+Gplusprime)/2);
            rhomarch(i,j)=U(1,1);
            umarch(i,j)=U(1,2)./rhomarch(i,j);
            vmarch(i,j)=U(1,3)./rhomarch(i,j);
            emarch(i,j)=U(1,4)./rhomarch(i,j)-umarch(i,j).^2/2-vmarch(i,j).^2/2;
            Pmarch(i,j)=emarch(i,j)*(gamma-1).*rhomarch(i,j); 
        end
        rhomarch(1,j)=rhomarch(2,j); rhomarch(M,j)=rhomarch(M-1,j);
        umarch(1,j)=umarch(2,j); umarch(M,j)=umarch(M-1,j);
        vmarch(1,j)=vmarch(2,j); vmarch(M,j)=vmarch(M-1,j);
        emarch(1,j)=emarch(2,j); emarch(M,j)=emarch(M-1,j);
        Pmarch(1,j)=Pmarch(2,j); Pmarch(M,j)=Pmarch(M-1,j);
        rho(:,j)=rhomarch(:,j);
        u(:,j)=umarch(:,j);
        v(:,j)=vmarch(:,j);
        P(:,j)=Pmarch(:,j);
        e(:,j)=emarch(:,j);
    end
    
    t=t+dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1:M;
j=1:M;
figure(1);
x=linspace(-(M-shift)*dx,(M-shift)*dx); 
y=linspace(-(M-shift)*dx,(M-shift)*dx);
image(x,y,rho(i,j),'CDataMapping','scaled');
colorbar
set(gca,'Yscale','lin','Ydir','normal');
set(gca, 'FontSize', 14);
title('Density')
xlabel('x')
ylabel('y')

figure(2); 
x=linspace(-(M-shift)*dx,(M-shift)*dx); 
y=linspace(-(M-shift)*dx,(M-shift)*dx);
image(x,y,u(i,j),'CDataMapping','scaled');
colorbar
set(gca,'Yscale','lin','Ydir','normal');
set(gca, 'FontSize', 14);
title('Velocity x-component')
xlabel('x')
ylabel('y')

figure(3); 
x=linspace(-(M-shift)*dx,(M-shift)*dx); 
y=linspace(-(M-shift)*dx,(M-shift)*dx);
image(x,y,v(i,j),'CDataMapping','scaled');
colorbar
set(gca,'Yscale','lin','Ydir','normal');
set(gca, 'FontSize', 14);
title('Velocity y-component')
xlabel('x')
ylabel('y')

figure(4); 
x=linspace(-(M-shift)*dx,(M-shift)*dx); 
y=linspace(-(M-shift)*dx,(M-shift)*dx);
image(x,y,P(i,j),'CDataMapping','scaled');
colorbar
set(gca,'Yscale','lin','Ydir','normal');
set(gca, 'FontSize', 14);
title('Pressure')
xlabel('x')
ylabel('y')

figure(5); 
x=linspace(-(M-shift)*dx,(M-shift)*dx); 
y=linspace(-(M-shift)*dx,(M-shift)*dx);
image(x,y,e(i,j),'CDataMapping','scaled');
colorbar
set(gca,'Yscale','lin','Ydir','normal');
set(gca, 'FontSize', 14);
title('Internal energy per unit mass')
xlabel('x')
ylabel('y')

function Fhllciplushalf = HLLCFlux(SL,SR,Sstar,FL,FR,FstarL,FstarR)
    if SL>=0
        Fhllciplushalf = FL;
    elseif SL<0 && Sstar>=0
        Fhllciplushalf = FstarL;
    elseif Sstar<0 && SR>=0
        Fhllciplushalf = FstarR;
    elseif SR<0
        Fhllciplushalf = FR;
    end
end

function [SL,SR,Sstar,FL,FR,FstarL,FstarR] = VectorsX(rhoL,rhoR,uL,uR,vL,vR,PL,PR)
        gamma=1.4;
        aL=(gamma*PL/rhoL)^0.5;
        aR=(gamma*PR/rhoR)^0.5;
        rhobar=(rhoL+rhoR)/2; abar=(aL+aR)/2;
        Ppvrs=(PL+PR)/2-(uR-uL)*rhobar*abar/2;
        P=max([0 Ppvrs]);
        if P<=PL
            qL=1;
        elseif P>PL
            qL=(1+(gamma+1)/(2*gamma)*(P/PL-1))^0.5;
        end
        if P<=PR
            qR=1;
        elseif P>PR
            qR=(1+(gamma+1)/(2*gamma)*(P/PR-1))^0.5;
        end
        SL=uL-aL*qL; SR=uR+aR*qR;
        Sstar=(PR-PL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL)-rhoR*(SR-uR));
        eL=PL/(gamma-1)/rhoL;
        EL=rhoL*(uL^2/2+eL);
        eR=PR/(gamma-1)/rhoR;
        ER=rhoR*(uR.^2/2+eR);
        UL=[rhoL rhoL*uL rhoL*vL EL];
        UR=[rhoR rhoR*uR rhoR*vR ER];
        FL=[rhoL*uL (rhoL*uL.^2+PL) rhoL*uL*vL uL*(EL+PL)];
        FR=[rhoR*uR (rhoR*uR.^2+PR) rhoR*uR*vR uR*(ER+PR)];
        UstarL=rhoL*(SL-uL)/(SL-Sstar)*[1 Sstar vL (EL/rhoL+(Sstar-uL)*(Sstar+PL/(rhoL*(SL-uL))))];
        UstarR=rhoR*(SR-uR)/(SR-Sstar)*[1 Sstar vR (ER/rhoR+(Sstar-uR)*(Sstar+PR/(rhoR*(SR-uR))))];
        FstarL=FL+SL*(UstarL-UL);
        FstarR=FR+SR*(UstarR-UR);
end

function [SL,SR,Sstar,GL,GR,GstarL,GstarR] = VectorsY(rhoL,rhoR,uL,uR,vL,vR,PL,PR)
        gamma=1.4;
        aL=(gamma*PL/rhoL)^0.5;
        aR=(gamma*PR/rhoR)^0.5;
        rhobar=(rhoL+rhoR)/2; abar=(aL+aR)/2;
        Ppvrs=(PL+PR)/2-(vR-vL)*rhobar*abar/2;
        P=max([0 Ppvrs]);
        if P<=PL
            qL=1;
        elseif P>PL
            qL=(1+(gamma+1)/(2*gamma)*(P/PL-1))^0.5;
        end
        if P<=PR
            qR=1;
        elseif P>PR
            qR=(1+(gamma+1)/(2*gamma)*(P/PR-1))^0.5;
        end
        SL=vL-aL*qL; SR=vR+aR*qR;
        Sstar=(PR-PL+rhoL*vL*(SL-vL)-rhoR*vR*(SR-vR))/(rhoL*(SL-vL)-rhoR*(SR-vR));
        eL=PL/(gamma-1)/rhoL;
        EL=rhoL*(vL^2/2+eL);
        eR=PR/(gamma-1)/rhoR;
        ER=rhoR*(vR.^2/2+eR);
        UL=[rhoL rhoL*uL rhoL*vL EL];
        UR=[rhoR rhoR*uR rhoR*vR ER];
        GL=[rhoL*vL rhoL*uL*vL (rhoL*vL.^2+PL) vL*(EL+PL)];
        GR=[rhoR*vR rhoR*uR*vR (rhoR*vR.^2+PR) vR*(ER+PR)];
        UstarL=rhoL*(SL-vL)/(SL-Sstar)*[1 uL Sstar (EL/rhoL+(Sstar-vL)*(Sstar+PL/(rhoL*(SL-vL))))];
        UstarR=rhoR*(SR-vR)/(SR-Sstar)*[1 uR Sstar (ER/rhoR+(Sstar-vR)*(Sstar+PR/(rhoR*(SR-vR))))];
        GstarL=GL+SL*(UstarL-UL);
        GstarR=GR+SR*(UstarR-UR);
end