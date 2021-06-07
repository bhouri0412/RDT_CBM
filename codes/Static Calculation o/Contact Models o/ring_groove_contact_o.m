function [fcg,mcg]=...
    ring_groove_contact_o(gcl,dz_g,zg,zr,alr,beta,beta_th,y1,y2,sigmag,z,Pk,omega,isu)

Npe=length(zr);
Nby=50;

fcg=zeros(1,Npe);
mcg=zeros(1,Npe);

h0=gcl/2+isu*(zg+dz_g-zr);
al=(alr-beta-beta_th);
al=al.*(abs(al)>1e-8)+1e-8*(abs(al)<1e-8);

minh=(h0-isu*y2*al).*(isu*al>0) ...
    +(h0-isu*y1*al).*(isu*al<0);
contact=minh<(4*sigmag);
%contact region
ncontact=nnz(contact);
if ncontact
    yc=(4*sigmag-h0(contact))./(-isu*al(contact));
    yc=min(yc,y2);
    yc=max(yc,y1);
    ya=y1*(isu*al(contact)<0)+yc.*(isu*al(contact)>0);
    ya=((y2-ya)>1e-8).*ya+(y2-1e-8).*((y2-ya)<1e-8);  %ya cannot excess y2
    yb=y2*(isu*al(contact)>0)+yc.*(isu*al(contact)<0);
    yb=((yb-y1)>1e-8).*yb+(y1+1e-8).*((yb-y1)<1e-8);    %yb cannot smaller than y1
    
    dy=(yb-ya)/Nby;
    dy2=repmat(dy,Nby,1);
    yy=repmat(ya,Nby+1,1)+repmat([0:Nby]',1,ncontact).*repmat(dy,Nby+1,1);
    h0_2=repmat(h0(contact),Nby+1,1);
    al_2=repmat(al(contact),Nby+1,1);
    h=h0_2-isu*al_2.*yy;
    %h=h.*(h>1e-8)+1e-8*(h<1e-8);
    Pcg=Pk*(omega-h/sigmag).^z.*(h/sigmag<=omega);
    dmcg=Pcg.*yy;
    fcg(contact)=sum(-isu*(Pcg(1:end-1,:)+Pcg(2:end,:))/2.*dy2); %in N/m
    mcg(contact)=sum(-isu*(dmcg(1:end-1,:)+dmcg(2:end,:))/2.*dy2);
end
        

end



