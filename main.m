clc,clear
%% Initialization, calculate Rw & Rl in every CC
N=5;  % CC number:5
Rw_min=2;  % wifi minimun average STA thr:2Mbps
W_num=[5,5,4,4,5];  % wifi STA num in different CCs
L_num=[10,10,10,10,10];  % initial LTE UE num in different CCs
W_num1=W_num;W_num2=W_num;  
L_num1=L_num;L_num2=L_num;  % other CC selection methods
Rw_avg=zeros(1,N);  % average STA thr in CC 1-N
Rl_avg=zeros(1,N);  % average user thr in CC 1-N
Rw_avg1=Rw_avg;Rw_avg2=Rw_avg;
Rl_avg1=Rl_avg;Rl_avg2=Rl_avg;  % other CC selection methods
W_temp=zeros(1,N);L_temp=zeros(1,N);  % temp value
% loc_x=[-20,25,-30,35,-40];
% loc_y=[-20,25,-30,35,-40];  % wifi AP location, LTE thr is also calculated by it
% Rw=zeros(1,N);  % wifi thr in CC 1-N
% Rl=zeros(1,N);  % LTE thr in CC 1-N
% for i=1:N
%     Rw(i)=wifinet(W_num(i),loc_x(i),loc_y(i));
%     Rl(i)=LTEnet(loc_x(i),loc_y(i));
% end  % calculate wifi & LTE thr in CC 1-N
load('Rw.mat');
load('Rl.mat');

%% user arrives, assign users in CC
lamda=2;  % user arrive at Possion distribution:lamda=2
Tmax=80;  % maximun simulation time:50s 
i=1;a=random('exponential',lamda); 
T(1)=round(a*10)/10; 
w(1)=T(1);  % initialization
L_min(1)=1;
L_min1(1)=1;
L_min2(1)=1;  % minimum LTE-U average user thr
L_chs(1)=1;
L_chs1(1)=1;
L_chs2(1)=1;  % the choosed LTE-U average user thr
W_min(1)=1;
W_min1(1)=1;
W_min2(1)=1;  % the choosed wifi average user thr
while(w(i)<Tmax)       
    for j=1:N
        W_temp(j)=1/(W_num(j)+L_num(j)+1)*Rw(j);
        L_temp(j)=1/(W_num(j)+L_num(j)+1)*Rl(j);  % calculate average thr in CC 1-N
        if W_temp(j)<Rw_min
            L_temp(j)=0;
        end
    end
    [ma,k]=max(L_temp);  
    L_num(k)=L_num(k)+1;  % select the maximum CC      
    k1=floor(rand*N)+1;
    L_num1(k1)=L_num1(k1)+1;  % random selction
    k2=mod(i,N)+1;
    L_num2(k2)=L_num2(k2)+1;  % circular selction
    for j=1:N
        Rw_avg(j)=Rw(j)/(W_num(j)+L_num(j));
        Rw_avg1(j)=Rw(j)/(W_num1(j)+L_num1(j));
        Rw_avg2(j)=Rw(j)/(W_num2(j)+L_num2(j));
        Rl_avg(j)=Rl(j)/(W_num(j)+L_num(j));
        Rl_avg1(j)=Rl(j)/(W_num1(j)+L_num1(j));
        Rl_avg2(j)=Rl(j)/(W_num2(j)+L_num2(j));
    end
    L_min(i)=min(Rl_avg);
    L_min1(i)=min(Rl_avg1);
    L_min2(i)=min(Rl_avg2); 
    L_chs(i)=Rl_avg(k);
    L_chs1(i)=Rl_avg1(k1);
    L_chs2(i)=Rl_avg2(k2);
    W_min(i)=min(Rw_avg);
    W_min1(i)=min(Rw_avg1);
    W_min2(i)=min(Rw_avg2);
    T(i)=random('exponential',lamda);  % user arriving interval follows Possion 
    T(i)=round(T(i)*10)/10; 
    w(i+1)=w(i)+T(i); i=i+1; % user arriving moment
end

%% plot the results
figure (1)
avgl=zeros(N,3);
for i=1:N
    avgl(i,1)=Rl_avg(i);
    avgl(i,2)=Rl_avg1(i);
    avgl(i,3)=Rl_avg2(i);
end
fx =bar(avgl);
ch = get(fx,'children');
set(fx(1), 'FaceColor',[1 0 0]);% Pink
set(fx(2), 'FaceColor',[0 1 0]); % Green
set(fx(3), 'FaceColor',[0 0 1]); % Green
set(gca,'XTickLabel',{'1','2','3','4','5'}) %设置x轴所代表大时间
xlabel('CC Index'),ylabel('Average User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('proposed','CC-R','CC-C')  %区分一下蓝色和红色分别代表什么
% title('LTE-U UE rate in every CC');      
figure (2)
avgw=zeros(N,3);
for i=1:N
    avgw(i,1)=Rw_avg(i);
    avgw(i,2)=Rw_avg1(i);
    avgw(i,3)=Rw_avg2(i);
end
fx=bar(avgw);
set(gca,'XTickLabel',{'1','2','3','4','5'}) %设置x轴所代表大时间
xlabel('CC Index'),ylabel('Average STA Throughput (Mbps)')  %设置x轴和y轴的名称
legend('proposed','CC-R','CC-C')  %区分一下蓝色和红色分别代表什么
ch = get(fx,'children');
set(fx(1), 'FaceColor',[1 0 0]);% Pink
set(fx(2), 'FaceColor',[0 1 0]); % Green
set(fx(3), 'FaceColor',[0 0 1]); % Green
hold on
x=get(gca,'xlim');
y=2;
plot(x,[y y],'k')
% title('WiFi STA rate in every CC');  
% x=1:length(L_min);
figure (3)
plot(L_min,'-or')
hold on
plot(L_min1,'-+g');
hold on
plot(L_min2,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Minimum LTE-U User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么
set(gca,'YLim',[4 8]);  %Y轴的数据显示范围
figure (4)
plot(W_min,'-or')
hold on
plot(W_min1,'-+g');
hold on
plot(W_min2,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Minimum WiFi User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么
x=get(gca,'xlim');
y=2;
plot(x,[y y],'k')  % threshold
set(gca,'YLim',[1.5 3]);  %Y轴的数据显示范围
figure (5)
plot(L_chs,'-or')
hold on
plot(L_chs1,'-+g');
hold on
plot(L_chs2,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Newly Arrived User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么

