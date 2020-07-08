clc,clear
%% Initialization, calculate Rw & Rl in every CC
N=5;  % CC number:5
Rw_min=2;  % wifi minimun average STA thr:1Mbps
W_num=[5,5,4,4,5];  % wifi STA num in different CCs
L_num=[10,10,10,10,10];  % initial LTE UE num in different CCs
W_num1=W_num;W_num2=W_num;  
L_num1=L_num;L_num2=L_num;  % other CC selection methods
Rl_avg=zeros(1,N);  % user average thr in CC 1-N
Rw_avg=zeros(1,N);  % user average thr in CC 1-N
Rl_rcv=zeros(1,N);  % user received thr in CC 1-N
Rw_rcv=zeros(1,N);  % user received thr in CC 1-N
Rw_avg1=Rw_avg;Rw_avg2=Rw_avg;
Rl_avg1=Rl_avg;Rl_avg2=Rl_avg;  % other CC selection methods
W_temp=zeros(1,N);L_temp=zeros(1,N);  % temp value
W_temp1=zeros(1,N);L_temp1=zeros(1,N);  % temp value
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

%% consider user type
UE_time=[2 1.2 0.8 0.6 0.4];  % set four types of users, occupying different channel time
UE_Uti=[3 1.5 1 0.6 0.2];  % User utility
Ltime=zeros(1,N);  % LTE user channel occupation time
for i=1:N
    for j=1:L_num(i)
        UE_type=floor(rand*length(UE_time))+1;  % UE type, random selected
        Ltime(i)=Ltime(i)+UE_time(UE_type);
    end
end

%% user arrives, assign users in CC
lamda=2;  % user arrive at Possion distribution:lamda=2
Tmax=60;  % maximun simulation time:50s 
% max_time=1;  % times to run
i=1;a=random('exponential',lamda); 
T(1)=round(a*10)/10; 
w(1)=T(1);  % initialization
L_Uti(1)=0;
L_Uti1(1)=0;
L_Uti2(1)=0;
L_Uti3(1)=0;
while(w(i)<Tmax)       
    for j=1:N
        UE_type=floor(rand*length(UE_time))+1;  % UE type, random selected
        alpha=UE_time(UE_type);  % random selected channel time
        W_temp(j)=1/(W_num(j)+Ltime(j)+alpha)*Rw(j);
        L_temp(j)=alpha/(W_num(j)+Ltime(j)+alpha)*Rl(j);  % calculate received thr in CC 1-N
        if W_temp(j)<Rw_min
            L_temp(j)=0;
        end
        W_temp1(j)=1/(W_num(j)+L_num(j)+1)*Rw(j);
        L_temp1(j)=1/(W_num(j)+L_num(j)+1)*Rl(j);  % calculate average thr in CC 1-N
        if W_temp1(j)<Rw_min
            L_temp1(j)=0;
        end
    end
    [ma,k]=max(L_temp);    
    Ltime(k)=Ltime(k)+alpha;  % select the maximum CC
    L_Uti(i)=ma*UE_Uti(UE_type);
    [ma1,k1]=max(L_temp1); 
    L_num(k1)=L_num(k1)+1; 
    L_Uti1(i)=ma1*UE_Uti(UE_type);
    k2=floor(rand*N)+1;
    L_num1(k2)=L_num1(k2)+1;  % random selction
    th1=Rl(k2)/(L_num1(k2)+W_num1(k2));
    L_Uti2(i)=th1*UE_Uti(UE_type);
    k3=mod(i,N)+1;
    L_num2(k3)=L_num2(k3)+1;  % circular selction
    th2=Rl(k3)/(L_num1(k3)+W_num1(k3));
    L_Uti3(i)=th2*UE_Uti(UE_type);   
    T(i)=random('exponential',lamda);  % user arriving interval follows Possion 
    T(i)=round(T(i)*10)/10; 
    w(i+1)=w(i)+T(i); i=i+1; % user arriving moment  
end
figure
x=1:length(L_Uti);
plot(x,L_Uti,'-^m');
hold on
plot(x,L_Uti1,'-or');
plot(x,L_Uti2,'-+g');
plot(x,L_Uti3,'-sb');
% grid on
xlabel('UE arriving number'),ylabel('User Utility')  %设置x轴和y轴的名称
legend('proposed-A','proposed','CC-R','CC-C')  %区分一下蓝色和红色分别代表什么