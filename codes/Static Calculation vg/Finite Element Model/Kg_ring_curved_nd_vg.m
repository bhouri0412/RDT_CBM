function [Kg]=Kg_ring_curved_nd_vg(Er,Gr,Iyr,Jtr,Rr,Telt,Nbe)
%Generate the element stiffness matrix using 3D curved beam element by
%integration (trapeze method)
%Initial geometry is assumed to be circular

%element length
Le=Rr*Telt;

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

%========= for twist displacements ===========%
Nt1=1-3*x.^2+2*x.^3;
Nt2=x-2*x.^2+x.^3;
Nt3=3*x.^2-2*x.^3;
Nt4=-x.^2+x.^3;

dNt1=-6*x+6*x.^2;
dNt2=1-4*x+3*x.^2;
dNt3=6*x-6*x.^2;
dNt4=-2*x+3*x.^2;

Nt=[Nt1;Nt2;Nt3;Nt4];
dNt=[dNt1;dNt2;dNt3;dNt4];
%=========================================%
Kge=zeros(16,16);
Kye=zeros(6,6);
Kzae=zeros(10,10);

for i=1:6
    for j=1:6
        dKij=(Telt^2*N(i,:)+d2N(i,:)).*(Telt^2*N(j,:)+d2N(j,:));
        Kye(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
    end
end

k=[1,2,3,1,2,4,5,6,3,4];

for i=1:10
    for j=1:10
        %i belongs to kz, j belongs to kz
        if mod(i,5)>=1&&mod(i,5)<=3&&mod(j,5)>=1&&mod(j,5)<=3
            dKij=d2N(k(i),:).*d2N(k(j),:)+Gr*Jtr*Telt^2/Er/Iyr*dN(k(i),:).*dN(k(j),:);
            Kzae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
        %i belongs to kz, j belongs to ka            
        elseif mod(i,5)>=1&&mod(i,5)<=3&&(mod(j,5)==0||mod(j,5)==4)
            dKij=Telt^2*d2N(k(i),:).*Nt(k(j),:)-Gr*Jtr*Telt^2/Er/Iyr*dN(k(i),:).*dNt(k(j),:);
            Kzae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
        %i belongs to ka, j belongs to kz
        elseif (mod(i,5)==0||mod(i,5)==4)&&mod(j,5)>=1&&mod(j,5)<=3
            dKij=Er*Iyr/Gr/Jtr*Nt(k(i),:).*d2N(k(j),:)-dNt(k(i),:).*dN(k(j),:);
            Kzae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
        %i belongs to ka, j belongs to ka
        else
            dKij=Er*Iyr/Gr/Jtr*Telt^2*Nt(k(i),:).*Nt(k(j),:)+dNt(k(i),:).*dNt(k(j),:);
            Kzae(i,j)=sum((dKij(:,1:end-1)+dKij(:,2:end))/2*dx);
        end     
        
    end
end

Kge([1:3,9:11],[1:3,9:11])=Kye;
Kge([4:8,12:16],[4:8,12:16])=Kzae;

Kg=zeros(8*(Nbe+1));

for i=1:Nbe
   %Filling the global stiffness matrix
   indg=(8*(i-1)+1):(8*(i+1));
   Kg(indg,indg)=Kg(indg,indg)+Kge;
end



