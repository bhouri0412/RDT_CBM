function [Ccs]=Cs_3rd_cont_vg(Tnod)

Nbnod=length(Tnod);

%Calculate the matrix form of the equation of continuity of the 3rd order
%derivative equations of contimuity were derived from the 3rd derivatives
%of the shape functions calculated in "Shape function 5th" Mathematica
%notebook

Ad2y=zeros(Nbnod-2); %initialization
Telt1=Tnod(2)-Tnod(1);
Telt2=Tnod(3)-Tnod(2);
Ad2y(1,1:2)=[9*(Telt1^2+Telt2^2),-3*Telt2^2]; %first line
Teltnm1=Tnod(Nbnod-1)-Tnod(Nbnod-2);
Teltn=Tnod(Nbnod)-Tnod(Nbnod-1);
Ad2y(Nbnod-2,Nbnod-3:Nbnod-2)=[-3*Teltnm1^2,9*(Teltnm1^2+Teltn^2)]; %last line
%filling the matrix
for i=2:Nbnod-3
    Teltim1=Tnod(i+1)-Tnod(i);
    Telti=Tnod(i+2)-Tnod(i+1);
    Ad2y(i,i-1:i+1)=[-3*Teltim1^2,9*(Teltim1^2+Telti^2),-3*Telti^2];
end

%Matrix defining the right hand side of Ad2y*d2y=B
By=zeros(Nbnod-2,2*Nbnod+2); %initialization
%filling the matrix
for i=2:Nbnod-1
    Teltim1=Tnod(i)-Tnod(i-1);
    Telti=Tnod(i+1)-Tnod(i);
    By(i-1,2*(i-1):2*i+3)=...
    [60,24*Teltim1,-120,36*(Teltim1-Telti),60,-24*Telti];
end

By(1,1)=3*Telt1^2; %completing first line
By(Nbnod-2,2*Nbnod+2)=3*Teltn^2; %completing last line

%Calculating second derivatives with the continuity constraint
Ccs=Ad2y\By;
end
