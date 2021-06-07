function [ybn]=bore_dist_static_o(inputs)
%computing bore distortion at ring nodes and current bore height

%Retrieving parmaters from input structure
Ab=inputs.bore.Ab; %bore distortiion magnitude
Phib=inputs.bore.Phib; %bore distortion phase
Tnod=inputs.FEM.Tnod; %angular position of nodes
Nbnod=inputs.FEM.Nbnod; %number of nodes
Tg=inputs.ring.Tg;
Nbe=inputs.FEM.Nbe;

uref=inputs.FEM.uref;
Telt=inputs.FEM.Telt;

ndb=length(Ab); %number of distortion orders

%transforming into column vector
Tnod=Tnod';
Tnodd=(0:2*pi/Nbe:2*pi)'+Tg;

%indexes for filling nodal displacement matrices
id1=1:3:3*Nbnod-2;
id2=2:3:3*Nbnod-1;
id3=3:3:3*Nbnod;

%initializing
ybn=zeros(3*Nbnod,1);

%bore distortion at node location
for i=2:ndb
    ybn(id1)=ybn(id1)+Ab(i)*sin(i*(Tnodd+Phib(i-1)))/uref;
    ybn(id2)=ybn(id2)+i*Ab(i)*cos(i*(Tnodd+Phib(i-1)))/uref*Telt;
    ybn(id3)=ybn(id3)-i^2*Ab(i)*sin(i*(Tnodd+Phib(i-1)))/uref*Telt^2;
end
%adding 0th order distortion
ybn(id1)=ybn(id1)+Ab(1)/uref;

end

