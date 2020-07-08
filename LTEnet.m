function [ Rl] = LTEnet( x1,y1 )
%输入为载波i用户所在的平均位置
%输出为LTE网络的用户平均吞吐量
%本函数根据信道容量公式，求出LTE网络的吞吐量
P=20;  % LTE-U SBS tx power:20dBm
g=10;  % channel gain:10
sigma=-95;  % noise power:-95dBm
B=20;  % carrier bandwidth:20MHz
x0=0;y0=0;  % location of SBS
d=sqrt((x0-x1)^2+(y0-y1)^2);  % distance between SBS and cluster center
alpha=5;  % pathloss coefficient in unlicensed band
pl=15.3+alpha*10*log10(d);  % path loss
Px=10^((P-pl)/10);  % received power in mW
N0=10^(sigma/10);  % noise power in mW
Rl=B*log2(1+g*Px/N0);  % data rate
% avg_thr=Rl/num;  % average user throughput
end

