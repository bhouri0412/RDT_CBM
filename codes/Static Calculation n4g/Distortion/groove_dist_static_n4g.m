function [zgnu,zgnl]=groove_dist_static_n4g(inputs)
%computing bore distortion at ring nodes and current bore height

%Retrieving parmaters from input structure
Agu=inputs.groove.Agu; %upper groove distortion magnitude
Phigu=inputs.groove.Phigu; %upper groove distortion phase
Agl=inputs.groove.Agl; %lower groove distortion magnitude
Phigl=inputs.groove.Phigl; %lower groove distortion phase
Tnod=inputs.FEM.Tnod; %angular position of nodes
Nbnod=inputs.FEM.Nbnod; %number of nodes
href=inputs.FEM.href;
Telt=inputs.FEM.Telt;
Tg=inputs.ring.Tg;
Nbe=inputs.FEM.Nbe;

ngdu=length(Agu); %number of distortion orders - upper flank
ngdl=length(Agl); %number of distortion orders - lower flank

%transforming into column vector
Tnod=Tnod';
Tnodd=(0:2*pi/Nbe:2*pi)'+Tg;

%indexes for filling nodal displacement matrices
id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

%initializing
zgnu=zeros(3*Nbnod,1);
zgnl=zeros(3*Nbnod,1);

%groove distortion at node location
%upper flank
for i=1:ngdu
    zgnu(id1)=zgnu(id1)+Agu(i)*sin((i+1)*(Tnodd+Phigu(i)))/href;
    zgnu(id2)=zgnu(id2)+(i+1)*Agu(i)*cos((i+1)*(Tnodd+Phigu(i)))/href*Telt;
    zgnu(id3)=zgnu(id3)-(i+1)^2*Agu(i)*sin((i+1)*(Tnodd+Phigu(i)))/href*Telt^2;
end

%lower flank
for i=1:ngdl
    zgnl(id1)=zgnl(id1)+Agl(i)*sin((i+1)*(Tnodd+Phigl(i)))/href;
    zgnl(id2)=zgnl(id2)+(i+1)*Agl(i)*cos((i+1)*(Tnodd+Phigl(i)))/href*Telt;
    zgnl(id3)=zgnl(id3)-(i+1)^2*Agl(i)*sin((i+1)*(Tnodd+Phigl(i)))/href*Telt^2;
end
end

