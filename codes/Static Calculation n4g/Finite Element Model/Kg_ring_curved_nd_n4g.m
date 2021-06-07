function [Kg]=Kg_ring_curved_nd_n4g(Er,Gr,Iyr,Jtr,Rr,Telt,Nbe,arm,ali,alo,isi,fref,href)
%Generate the element stiffness matrix using 3D curved beam element by
%integration (trapeze method)
%Initial geometry is assumed to be circular

%element length
Le=Rr*Telt;

aleff=ali*(isi==1)-alo*(isi==0);

%discretization inside the element
Npe=1000;
dx=1/Npe;
x=0:dx:1;

%Shape functions and their derivatives
%========= for y and z displacements ===========%
N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=x-6*x.^3+8*x.^4-3*x.^5;
N3=(x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2;
N4=10*x.^3-15*x.^4+6*x.^5;
N5=-4*x.^3+7*x.^4-3*x.^5;
N6=(x.^3)/2-x.^4+x.^5/2;

dN1=-30*x.^2+60*x.^3-30*x.^4;
dN2=1-18*x.^2+32*x.^3-15*x.^4;
dN3=x-9*(x.^2)/2+6*x.^3-5*(x.^4)/2;
dN4=30*x.^2-60*x.^3+30*x.^4;
dN5=-12*x.^2+28*x.^3-15*x.^4;
dN6=3*(x.^2)/2-4*x.^3+5*x.^4/2;

d2N1=-60*x+180*x.^2-120*x.^3;
d2N2=-36*x+96*x.^2-60*x.^3;
d2N3=1-9*x+18*x.^2-10*x.^3;
d2N4=60*x-180*x.^2+120*x.^3;
d2N5=-24*x+84*x.^2-60*x.^3;
d2N6=3*x-12*x.^2+10*x.^3;

N=[N1;N2;N3;N4;N5;N6];
dN=[dN1;dN2;dN3;dN4;dN5;dN6];
d2N=[d2N1;d2N2;d2N3;d2N4;d2N5;d2N6];

%=========================================%

Kge=zeros(24,24);

Kye=zeros(6,6);
Kzae=zeros(6,6);
Kaze=zeros(6,6);
Kae=zeros(6,6);
Kze=zeros(6,6);

Kzf=zeros(6,6);

for i=1:6
    for j=1:6
        dKij=(Telt^2*N(i,:)+d2N(i,:)).*(Telt^2*N(j,:)+d2N(j,:));
        Kye(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

for i=1:6
    for j=1:6
        dKij=d2N(i,:).*d2N(j,:)+Gr*Jtr*Telt^2/Er/Iyr*dN(i,:).*dN(j,:);
        Kze(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

for i=1:6
    for j=1:6
        dKij=Telt^2*d2N(i,:).*N(j,:)-Gr*Jtr*Telt^2/Er/Iyr*dN(i,:).*dN(j,:);
        Kzae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

for i=1:6
    for j=1:6
        dKij=Er*Iyr/Gr/Jtr*N(i,:).*d2N(j,:)-dN(i,:).*dN(j,:);
        Kaze(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

for i=1:6
    for j=1:6
        dKij=Er*Iyr/Gr/Jtr*Telt^2*N(i,:).*N(j,:)+dN(i,:).*dN(j,:);
        Kae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

for i=1:6
   for j=1:6
    dKij=N(i,:).*N(j,:);
    Kzf(i,j)=-Le^4*fref/Er/Iyr/href*sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
   end  
end

Kge([1:3,13:15],[1:3,13:15])=Kye;
Kge([4:6,16:18],[4:6,16:18])=Kze;
Kge([4:6,16:18],[7:9,19:21])=Kzae;
Kge([7:9,19:21],[4:6,16:18])=Kaze;
Kge([7:9,19:21],[7:9,19:21])=Kae;
Kge([4:6,16:18],[10:12,22:24])=Kzf;

Kge(10,4)=1;
Kge(10,7)=-aleff/Rr;
Kge(11,5)=1;
Kge(11,8)=-aleff/Rr;
Kge(12,6)=1;
Kge(12,9)=-aleff/Rr;

Kge(22,16)=1;
Kge(22,19)=-aleff/Rr;
Kge(23,17)=1;
Kge(23,20)=-aleff/Rr;
Kge(24,18)=1;
Kge(24,21)=-aleff/Rr;

Kg=zeros(12*(Nbe+1));

for i=1:Nbe
   %Filling the global stiffness matrix
   indg=(12*(i-1)+1):(12*(i+1));
   Kg(indg,indg)=Kg(indg,indg)+Kge;
end



