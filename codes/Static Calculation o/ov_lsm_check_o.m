function [U]=ov_lsm_check_o(ov,Nbe,Npe)
Tri=ov(:,1);
rfsim=ov(:,2);
Nbnod=Nbe+1;
%%5th interpolation of input data using LS method
Inod=1:Npe:Npe*Nbe+1;
dTrfs=(Tri(end)-Tri(1))/(Npe*Nbe);
Trfs=(Tri(1):dTrfs:Tri(length(Tri)))';
Teltfs=(Trfs(end)-Trfs(1))/Nbe;
Tnodfs=Trfs(Inod);

rfsi=interp1(Tri,rfsim,Trfs); %linear reinterpolation
[Ccs]=Cs_3rd_cont_o(Tnodfs);
[U]=ls_rfsi_o(rfsi,Nbnod,Teltfs,Ccs);

%figure(1)
%plot(Tri,(rfsim+Rr)*10^3,'Color',[0 0.8 0],'LineWidth',1.5);
%x=0:0.001:1;

%for k=1:Nbe
%    ovkd=U(3*k-2)*(1-10*x.^3+15*x.^4-6*x.^5)+U(3*k-1)*(Teltfs*(x-6*x.^3+8*x.^4-3*x.^5))+U(3*k)*(Teltfs^2*((x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2))+...
%        U(3*k+1)*(10*x.^3-15*x.^4+6*x.^5)+U(3*k+2)*(Teltfs*(-4*x.^3+7*x.^4-3*x.^5))+U(3*k+3)*(Teltfs^2*((x.^3)/2-x.^4+x.^5/2));
%    ovkd=(Rr+ovkd)*10^3;
%    y=Teltfs*x+(k-1)*Teltfs+Tri(1);
%    hold on
%    plot(y,ovkd);
%    min(ovkd)
%    grid on
%end
end