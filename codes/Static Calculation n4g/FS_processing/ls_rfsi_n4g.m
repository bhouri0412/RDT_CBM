function [rnod]=ls_rfsi_n4g(rfsi,Nbnod,Telt,Ccs)
%interpolation of input free shape radius with 5th order polynomial 
%spline. The spline is fitted on the input data using a least square method

Nr=length(rfsi)-1;
Nbe=Nbnod-1;
Npe=Nr/Nbe;

%======= Full matrix before consideration of continuity constraints ======%

%discretization inside the element
dx=1/Npe;
x=(0:dx:1)';

%Shape functions
N1=1-10*x.^3+15*x.^4-6*x.^5;
N2=Telt*(x-6*x.^3+8*x.^4-3*x.^5);
N3=Telt^2*((x.^2)/2-3*(x.^3)/2+3*(x.^4)/2-(x.^5)/2);
N4=10*x.^3-15*x.^4+6*x.^5;
N5=Telt*(-4*x.^3+7*x.^4-3*x.^5);
N6=Telt^2*((x.^3)/2-x.^4+x.^5/2);

%initialization of matrix used in least square method
Xlsf=zeros(Nr+1,3*Nbnod);
%derivative of approximating function
for i=2:Nbnod-1 %iterations on the nodes
    int=1+(i-2)*Npe:1+i*Npe;
    for j=1:3 %iteration on the nodal DOFs
        jDoFs=j+3*(i-1);
        switch j
            case 1
        Xlsf(int,jDoFs)=[N4;N1(2:Npe+1)];
            case 2
        Xlsf(int,jDoFs)=[N5;N2(2:Npe+1)];
            case 3
        Xlsf(int,jDoFs)=[N6;N3(2:Npe+1)];        
        end
    end
end
%completing remaining terms of the matrix for node 1 and n
int=1:1+Npe;
Xlsf(int,1:3)=[N1,N2,N3];
int=Nr+1-Npe:Nr+1;
Xlsf(int,3*Nbnod-2:3*Nbnod)=[N4,N5,N6];

%least square solution
Als=(Xlsf')*Xlsf;
Bls=(Xlsf')*rfsi;
rnod=Als\Bls;

%=========================================================================%
const=1;
if const
%=================== Adding the continuity constraints ===================%

%Retriving what is relevant to the independent parameters
Xls=zeros(Nr+1,2*Nbnod+2);
Xls(:,[1,2*Nbnod+2])=[Xlsf(:,3),Xlsf(:,3*Nbnod)]; 
Xls(:,2:2:2*Nbnod)=Xlsf(:,1:3:3*Nbnod-2);
Xls(:,3:2:2*Nbnod+1)=Xlsf(:,2:3:3*Nbnod-1);

%contrbution of constraint equations
for i=2:Nbnod-1
    id2rn=3*i;
    Xls=Xls+Xlsf(:,id2rn)*Ccs(i-1,:);
end

%=========================================================================%

%least square solution
Als=(Xls')*Xls;
Bls=(Xls')*rfsi;
rnodp=Als\Bls;

%completing displacements with continuity constraints
d2rnod=Ccs*rnodp;

rnod=zeros(3*Nbnod,1);
rnod(1:3:3*Nbnod-2)=rnodp(2:2:2*Nbnod);
rnod(2:3:3*Nbnod-1)=rnodp(3:2:2*Nbnod+1);
rnod(3:3:3*Nbnod)=[rnodp(1);d2rnod;rnodp(2*Nbnod+2)];
end
end

