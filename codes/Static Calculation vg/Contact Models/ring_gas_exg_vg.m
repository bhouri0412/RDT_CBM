function [frgas,fzgas,mgas]=ring_gas_exg_vg(Pi,Pu,Pd,arm,aui,lfu,ali,lfl,hui,hli,huo,hlo,rbn,rb1,rb2,a11,a12,a21,a22,alr)

z0_1=(rbn-(a11+alr)/2/a12);
z0_1=min(z0_1,rb2+rbn);
z0_1=max(z0_1,-rb1+rbn);
z0_2=(rbn-(a21+alr)/2/a22);
z0_2=min(z0_2,rb2+rbn);
z0_2=max(z0_2,-rb1+rbn);
checktwopoints=(z0_1<rbn).*(z0_2>rbn);

if nnz(checktwopoints)
    warning('bad ring profile');
    z0=rbn;
else
    z0=z0_1.*(z0_1<rbn).*(z0_2<rbn)+z0_2.*(z0_2>rbn).*(z0_1>rbn)+rbn*(z0_1>=rbn).*(z0_2<=rbn);
end

frgas=Pi.*(hui+hli)-Pu.*(huo-z0)-Pd.*(hlo+z0);
fzgas=-Pu*(arm+aui-lfu)+Pd*(arm+ali-lfl)-Pi*(ali-aui);
mgas=-1/2*Pi*hui^2+1/2*Pi*hli^2+1/2*Pu.*(huo^2-z0.^2)-1/2*Pd.*(hlo^2-z0.^2)...
    -1/2*Pu*(arm+lfu-aui)*(arm+aui-lfu)+1/2*Pd*(arm+lfl-ali)*(arm+ali-lfl)+1/2*Pi*(ali+aui)*(ali-aui);
end