clear all;
close all;

%Définition des constantes et initialisation

N=5;
To=100;
seuil=0.1;
h=0.001;
m=10;
g=9.8;
V=4;
C=0.008;
K=0.5*1.18*0.07*C;
Drag=K*V^2;
Bmax=0.4;
Tmax=80;
Lmax=787;
Vmax=15;
% a=0.45;
margeGamma=1/5;


nbEtapes=30 ;
path=zeros(3,nbEtapes);
numCible=1;


%Trajectoire de référence
Npoints=1000*nbEtapes;
GammaRef=zeros(1,Npoints);
PsiRef=zeros(1,Npoints);
GammaRef(1:Npoints/2)=0.4;
GammaRef(Npoints/2:Npoints)=0.2;

PsiRef(1:Npoints/2)=-pi/2;
PsiRef(Npoints/2:Npoints)=-pi/2+pi*(0:Npoints/2)/(Npoints/2);

VRef=zeros(3,Npoints);
MRef=zeros(3,Npoints);
MRef(:,1)=[0,0,0].';
compteur=0;
s=0;

for l = 1:(Npoints-1)
    VRef(:,l)=V*[cos(GammaRef(l))*sin(PsiRef(l)),cos(GammaRef(l))*cos(PsiRef(l)),sin(GammaRef(l))].';
    MRef(:,l+1)=MRef(:,l)+h*VRef(:,l);
    compteur=compteur+1;
    if(compteur==(Npoints/nbEtapes))
        compteur=0;
        s=s+1;
        path(:,s)=MRef(:,l+1);
    end
end
path(:,nbEtapes)=MRef(:,Npoints);

%Position initiale du drone
M=[0.1,0,0].';

%Lancement du drone

Gamma=0;
Psi=-pi/2;
VM=V*[cos(Gamma)*sin(Psi),cos(Gamma)*cos(Psi),sin(Gamma)].';

%Positions et vitesses relatives initiales
e=path(:,1)-M;
D=norm(e);

%Vecteurs de stockage des données

stop=To/h;
TrajM=zeros(3,stop);
Dist=zeros(1,stop);
GammaList=zeros(1,stop);
PsiList=zeros(1,stop);

Thrust=zeros(1,stop);
Lift=zeros(1,stop);
Bank=zeros(1,stop);
Dmax=0;

TrajM(:,1)=M;
Dist(1)=D;
GammaList(1)=Gamma;
PsiList(1)=Psi;

%Itération
%Méthode d'intégration utilisée: Euler explicite de pas h
for k = 1: stop
    e=path(:,numCible)-M;
    D=norm(e);
    
    % Détermination Dmax
    if (D>Dmax)
        Dmax=D;
    end
    if(D<seuil)
        if(numCible<nbEtapes)
            numCible=numCible+1;
            e=path(:,numCible)-M;
            D=norm(e);
        else
            kmax=k;
            break;
        end
    end
    
    %     e2=zeros(3,1);
    %     e3=zeros(3,1);
    %     e4=zeros(3,1);
    etot=e;
    %         if(numCible<nbEtapes)
    %             e2=path(:,numCible+1)-M;
    %             etot=etot+a*(e2.'*e)*e2/abs(e2.'*e);
    %             % etot=etot+a*e2;
    %         if(numCible<nbEtapes-1)
    %             e3=path(:,numCible+2)-M;
    %             etot=etot+a^2*(e3.'*e)*e3/abs(e3.'*e);
    %             % etot=etot+a^2*e3;
    %         end
    %         if(numCible<nbEtapes-2)
    %             e4=path(:,numCible+3)-M;
    %             etot=etot+a^3*(e4.'*e)*e4/abs(e4.'*e);
    %             % etot=etot+a^3*e4;
    %         end
    %     end
    Dtot=norm(etot);
    
    
    %Base tournante
    u=[cos(Gamma)*sin(Psi),cos(Gamma)*cos(Psi),sin(Gamma)];
    v=[cos(Psi),-sin(Psi),0];
    w=[sin(Gamma)*sin(Psi),sin(Gamma)*cos(Psi),-cos(Gamma)];
    
    dPsi=h*N*V*v*(etot)/(cos(Gamma)*Dtot^3); % oscillations T L Gamma
    dGamma=-h*N*V*w*(etot)/(Dtot^3); % faibles oscillations Bank
    
    %     dPsi=h*N*V*v*(etot)/(cos(Gamma)*Dtot^2); % pas oscillations T L Gamma
    %     dGamma=-h*N*V*w*(etot)/(Dtot^2); % oscillations Bank
    
    % Eviter la singularité Gamma=pi/2
    if(Gamma>=pi/2-margeGamma && Gamma<=pi/2)
        Gamma=pi/2-margeGamma;
    end
    if(Gamma<=pi/2+margeGamma && Gamma>pi/2)
        Gamma=pi/2+margeGamma;
    end
    
    Thrust(k)=Drag+m*g*sin(Gamma);
    Lc=m*g*cos(Gamma)+m*V*dGamma/h;
    Ls=m*V*cos(Gamma)*dPsi/h;
    Lift(k)=sqrt(Lc^2+Ls^2);
    Bank(k)=asin(Ls/Lift(k));
    sat=0;
    
    % saturations
    if (Thrust(k)>Tmax)
        Thrust(k)=Tmax;
        sat=sat+1;
    end
    
    if (Lift(k)>Lmax)
        Lift(k)=Lmax;
        Bank(k)=asin(Ls/Lift(k));
        sat=sat+1;
    end
    
    if (abs(Bank(k))>Bmax)
        Bank(k)=sign(Bank(k))*Bmax;
        sat=sat+1;
        
    end
    
    if(sat>=1)
        dPsi=h*Lift(k)*sin(Bank(k))/(m*V*cos(Gamma));
        dGamma=-h*g*cos(Gamma)/V+h*Lift(k)*cos(Bank(k))/(m*V);
        dV=h*(Thrust(k)-D)/m-h*g*sin(Gamma);
        V=V+dV;
        if (V>Vmax)
            V=Vmax;
        end
        Drag=K*C*V^2;
    end
    
    Psi=Psi+dPsi;
    Gamma=Gamma+dGamma;
    
    M=M+h*VM;
    VM=V*[cos(Gamma)*sin(Psi),cos(Gamma)*cos(Psi),sin(Gamma)].';
    
    GammaList(k+1)=Gamma;
    PsiList(k+1)=Psi;
    TrajM(:,k+1)=M;
    Dist(k+1)=D;
    
    if(k==stop)
        kmax=stop;
    end
end

%Visualisation

%Trajectoire
figure
plot3(TrajM(1,(1:kmax)),TrajM(2,(1:kmax)),TrajM(3,(1:kmax)),'r');
hold on
plot3(path(1,:),path(2,:),path(3,:),'+');
hold on
plot3(MRef(1,:),MRef(2,:),MRef(3,:));
grid on
title('Trajectoire drone (rouge) . Plan de vol (+). Trajectoire réf (bleu)');
xlabel('x');
ylabel('y');
zlabel('z');

figure
subplot(326)
%Distance cible-drone
% figure
plot(h*(1:kmax),Dist(1:kmax),'r');
title('Distance entre le drone et l étape courante.');
xlabel('Temps');
ylabel('D');
ylim([0 1.1*Dmax])

subplot(321)
%Pitch
% figure
plot(h*(1:kmax),GammaList(1:kmax),'r');
title('Gamma: réel(rouge) et réf (bleu)');
xlabel('Temps');
ylabel('Gamma');
hold on
plot((h*kmax/Npoints)*(1:Npoints),GammaRef);

subplot(323)
%Yaw
% figure
plot(h*(1:kmax),PsiList(1:kmax),'r');
title('Psi: réel(rouge) et réf (bleu)');
xlabel('Temps');
ylabel('Psi');
hold on
plot((h*kmax/Npoints)*(1:Npoints),PsiRef);

subplot(322)
%Thrust
% figure
plot(h*(1:kmax),Thrust(1:kmax));
title('Required Thrust');
xlabel('Temps');
ylabel('Thrust');
ylim([0 1.1*Tmax])

subplot(324)
%Lift
% figure
plot(h*(1:kmax),Lift(1:kmax));
title('Required Lift');
xlabel('Temps');
ylabel('Lift');
ylim([0 Lmax])

subplot(325)
%Bank
% figure
plot(h*(1:kmax),Bank(1:kmax));
title('Required Bank');
xlabel('Temps');
ylabel('Bank');
ylim([0 1.1*Bmax])






