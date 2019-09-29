function [ res ] = as_subcal2( A,B )
%AS_SUBCAL2 Summary of this function goes here
%   Detailed explanation goes here
[dm,dn]=size(A);
Q1=CS_gene(dm);
Q2=CS_gene(dn);
m=dm;
n=dn;
F1=dftmtx(m)/sqrt(m);
F2=dftmtx(n)/sqrt(n);
D1=F1'*Q1*F1;
D2=F2'*Q2*F2;
U1=F1;
U2=F2;
dA=diag(D1);
dB=diag(D2);
%% ’≈¡øº∆À„
FA=F1*A*F2;
FB=F1*B*F2;
C = FA./FB;
dA_u = (dA);
dB_u = dB;
angel_C = angle(C);
angel_A = dA_u*ones(1,dn);
angel_Au = angle(angel_A);
angel_B = dB_u*ones(1,dm);
angel_Bu = angle(angel_B);
angel_Bu = angel_Bu.';

% line = min(round(dm*0.25),round(dn*0.25));
% maskx = linspace(-dm/2,dm/2,dm)'*ones(1,dn);
% masky = ones(dm,1)*linspace(-dn/2,dn/2,dn);
% mask = (sqrt(maskx.^2+masky.^2)<line);
mask = (sqrt(angel_Au.^2+angel_Bu.^2)<0.5*pi);
angel_C = medfilt2(angel_C);
angel_x = angel_Au.*mask;
angel_y = angel_Bu.*mask;
angel_CC = angel_C.*mask;
%     angel_x = angel_Au(floor(dm/2)-line:floor(dm/2)+line,floor(dn/2)-line:floor(dn/2)+line);
%     angel_y = angel_Bu(floor(dm/2)-line:floor(dm/2)+line,floor(dn/2)-line:floor(dn/2)+line);
%     angel_CC = angel_C(floor(dm/2)-line:floor(dm/2)+line,floor(dn/2)-line:floor(dn/2)+line);
angel_x = reshape(angel_x,dm*dn,1);
angel_y = reshape(angel_y,dm*dn,1);
angel_CC = reshape(angel_CC,dm*dn,1);
angel_res = [angel_x,angel_y];
res = regress(angel_CC,angel_res);

end

