function [res] = ana10_ywt(img1,img2 ,ite_num)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
data1 = img1;
data2 = img2;
[dm,dn]=size(data1);
A = ifx(ify(data1));
B = ifx(ify(data2));
%%
C = A.*conj(B);
data3 = iifx(iify(C));
[mdata,ind1]=max(abs(data3));
[~,ind2] = max(mdata);
ind3 = ind1(ind2);
res1=[ind3-(dm+1)/2,ind2-(dn+1)/2+1];
res1=round(res1);
%% original data process
dx = res1(1);
dy = res1(2);
data3 = data1;
data4 = data2;
if(dx>0)
    data3 = data3(dx+1:end,:);
    data4 = data4(1:end-dx,:);
end
if(dx<0)
    data3 = data3(1:end+dx,:);
    data4 = data4(1-(dx):end,:);
end
if(dy>0)
    data3 = data3(:,dy+1:end);
    data4 = data4(:,1:end-dy);
end
if(dy<0)
    data3 = data3(:,1:end+dy);
    data4 = data4(:,1-dy:end);
end
[dm,dn]=size(data3);

m=dm;
n=dn;
% Q1(:,2:dm)=Q1(:,1:dm-1);
% Q1(dm,1)=1;
% Q1(1,1)=0;
% Q2(:,2:dn)=Q2(:,1:dn-1);
% Q2(dn,1)=1;
% Q2(1,1)=0;
Q1 = CS_gene(dm);
Q2 = CS_gene(dn);
F1=dftmtx(m)/sqrt(m);
F2=dftmtx(n)/sqrt(n);
D1=F1'*Q1*F1;
D2=F2'*Q2*F2;
dA=diag(D1);
dB=diag(D2);
A=data3;
B=data4;
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
% angel_C = medfilt2(angel_C);
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
res2 = regress(angel_CC,angel_res);
res=zeros(ite_num,2);
res(1,:)=res1+res2';
for i=2:ite_num
%   data4n=F1*(diag(dA.^res2(1)))*F1'*data4*F2*(diag(dB.^res2(2)))*F2';
    data4n=(Q1')^(res2(1))*data4*(Q2)^(res2(2));
    data3n=data3(2:end-1,2:end-1);
    data4n=data4n(2:end-1,2:end-1);
    resn=as_subcal2(data3n,data4n);
    res2=res2+resn;
    res(i,:)=res1+res2';
end


end

