clc,clear
%% To simulate and compare QoS & AC, UCCS, RS and CS algorithms
% Users are devided into two categories, QoS and BE users, with propotion 1:3
% QoS users have specific thr and delay requirements, we set it as 7,8,9,10 Mbps, devided with equal probability
% In fact, We did not consider the change of WiFi users. No Rw nor WiFi users num
% In UCCS, RS and CS, they did not have any different treatments of QoS users, but regard them as the same.
% In UCCS, RS and CS, they did not have any access control, thus the performance of WiFi users cannot be guaranteed.

%% Initialization, calculate Rw & Rl on every CC

N=5;  % CC number:5
Rw_min=2;  % wifi minimun average STA thr requirement:2Mbps
portion=0.4;  % the proportion of QoS users in all LTE-U users: 0.25
W_num=[5,5,4,4,5];  % initial wifi STA num in different CCs
B_num=[6,6,6,6,6];  % initial LTE BE UE num, initially all are BE users
Q_tim=[0,0,0,0,0];  % initial LTE QoS UE channel time for QoS & AC, initially no QoS users 
Q_num1=[0,0,0,0,0];  % initial LTE QoS UE num for UCCS, RS and CS, initially no QoS users 
Rq_req=[7,8,9,10];  % LTE QoS UE minimun thr requirements: Mbps
Del_req=[20,50,80,100];  % LTE QoS UE maximum delay requirements: ms
W_num1=W_num;W_num2=W_num; W_num3=W_num;
B_num1=B_num;B_num2=B_num; B_num3=B_num;
Q_num2=Q_num1;Q_num3=Q_num1;  % for UCCS, RS and CS CC selection methods

% loc_x=[-20,25,-30,35,-40];
% loc_y=[-20,25,-30,35,-40];  % wifi AP location, LTE thr is also calculated by it
% Rw=zeros(1,N);  % wifi thr in CC 1-N
% Rl=zeros(1,N);  % LTE thr in CC 1-N
% for i=1:N
%     Rw(i)=wifinet(W_num(i),loc_x(i),loc_y(i));
%     Rl(i)=LTEnet(loc_x(i),loc_y(i));
% end  
load('Rw.mat');
load('Rl.mat');  % calculate wifi & LTE thr in CC 1-N

Rw_avg=zeros(1,N);  
Rw_avg1=Rw_avg;Rw_avg2=Rw_avg;Rw_avg3=Rw_avg;  % average wifi STA thr on CC 1-N
Rb_avg=zeros(1,N);  
Rb_avg1=Rb_avg;Rb_avg2=Rb_avg;Rb_avg3=Rb_avg;  % average LTE BE UE thr on CC 1-N 
for j=1:N
    Rb_avg(j)=1/(W_num(j)+B_num(j))*(1-Q_tim(j))*Rl(j);
    Rb_avg1(j)=1/(W_num1(j)+B_num1(j))*Rl(j);
    Rb_avg2(j)=1/(W_num2(j)+B_num2(j))*Rl(j);
    Rb_avg3(j)=1/(W_num3(j)+B_num3(j))*Rl(j);
    Rw_avg(j)=1/(W_num(j)+B_num(j))*(1-Q_tim(j))*Rw(j);
    Rw_avg1(j)=1/(W_num1(j)+B_num1(j))*Rw(j);
    Rw_avg2(j)=1/(W_num2(j)+B_num2(j))*Rw(j);
    Rw_avg3(j)=1/(W_num3(j)+B_num3(j))*Rw(j);  % initialization
end
Rq_thr=zeros(1,N);
Rq_thr1=Rq_thr;Rq_thr2=Rq_thr;Rq_thr3=Rq_thr;  % LTE QoS UE thr on CC 1-N
W_temp=zeros(1,N);B_temp=zeros(1,N); Q_temp=zeros(1,N); % temp value

%% users arrive, assign users on CC

B_min(1)=0;
B_min1(1)=0;
B_min2(1)=0;
B_min3(1)=0;  % minimum LTE-U BE UE average user thr
B_chs(1)=0;
B_chs1(1)=0;
B_chs2(1)=0;
B_chs3(1)=0;  % the selected LTE-U BE UE average user thr
Q_chs(1)=0;
Q_chs1(1)=0;
Q_chs2(1)=0;
Q_chs3(1)=0;  % the selected LTE-U QoS UE average user thr
W_min(1)=0;
W_min1(1)=0;
W_min2(1)=0; 
W_min3(1)=0;  % the selected wifi average user thr

lamda=2;  % users arrive at Possion distribution:lamda=2
Tmax=150;  % maximun simulation time:150s 
Rq(1)=1;  % thr requirement of the arrived user
i=1;a=random('exponential',lamda);  
T(1)=round(a*10)/10; 
w(1)=T(1);  % initialization, arrival of the first LTE-U user

while(w(i)<Tmax)       
 % user type selection   
    temp=rand;  % rand number
    if temp < portion  % the arrived user is QoS
        m=floor(rand*4)+1;  
        Rq(i)=Rq_req(m);  % random select the thr requirement
    else
        Rq(i)=0;  % thr req of BE user
    end
  % for QoS & AC algorithm    
    if Rq(i) >0   % for QoS user
        for j=1:N
            Q_temp(j)=Q_tim(j)+Rq(i)/Rl(j);
            W_temp(j)=1/(W_num(j)+B_num(j))*(1-Q_temp(j))*Rw(j);
            B_temp(j)=1/(W_num(j)+B_num(j))*(1-Q_temp(j))*Rl(j);  % calculate average thr in CC 1-N
            if W_temp(j)<Rw_min
                B_temp(j)=0;
            end
        end
        [ma,k]=max(B_temp);  % find the CC that can provide the maximum thr for BE users
        if B_temp(k)==0
            Q_chs(i)=0;  % reject the access of such a QoS user
        else
            Q_tim(k)=Q_tim(k)+Rq(i)/Rl(k); % assign the QoS user on the selected CC
            Q_chs(i)=Rq(i);  
            Rb_avg(k)=1/(W_num(k)+B_num(k))*(1-Q_tim(k))*Rl(k);
            Rw_avg(k)=1/(W_num(k)+B_num(k))*(1-Q_tim(k))*Rw(k);  % update the change caused by the new QoS user  
        end
        B_chs(i)=0;
        B_min(i)=min(Rb_avg);
        W_min(i)=min(Rw_avg);  % the min thr of BE and WiFi users
    else   % for BE user
        for j=1:N
            W_temp(j)=1/(W_num(j)+B_num(j)+1)*(1-Q_tim(j))*Rw(j);
            B_temp(j)=1/(W_num(j)+B_num(j)+1)*(1-Q_tim(j))*Rl(j);  % calculate average thr in CC 1-N
            if W_temp(j)<Rw_min
                B_temp(j)=0;
            end
        end
        [ma,k]=max(B_temp);  % find the CC that can provide the maximum thr for BE users
        if B_temp(k)==0
            B_chs(i)=0;  % reject the access of such a BE user
        else
            B_num(k)=B_num(k)+1;  % assign the BE user on the selected CC
            B_chs(i)=Rb_avg(k);
            Rb_avg(k)=1/(W_num(k)+B_num(k))*(1-Q_tim(k))*Rl(k);
            Rw_avg(k)=1/(W_num(k)+B_num(k))*(1-Q_tim(k))*Rw(k);  % update the change caused by the new BE user           
        end
        Q_chs(i)=0;
        B_min(i)=min(Rb_avg);
        W_min(i)=min(Rw_avg);  % the min thr of BE and WiFi users
    end
  % for UCCS, RS and CS algorithms    
    for j=1:N
        W_temp(j)=1/(W_num1(j)+B_num1(j)+Q_num1(j)+1)*Rw(j);
        B_temp(j)=1/(W_num1(j)+B_num1(j)+Q_num1(j)+1)*Rl(j);  % calculate average thr in CC 1-N
        if W_temp(j)<Rw_min
            B_temp(j)=0;
        end
    end
    [ma1,k1]=max(B_temp);  % UCCS
    k2=floor(rand*N)+1;  % RS   
    k3=mod(i,N)+1;  % CS  
    if Rq(i)>0   % for QoS user
        Q_num1(k1)=Q_num1(k1)+1;
        Q_num2(k2)=Q_num2(k2)+1; 
        Q_num3(k3)=Q_num3(k3)+1; 
    else  % for BE user
        B_num1(k1)=B_num1(k1)+1;
        B_num2(k2)=B_num2(k2)+1; 
        B_num3(k3)=B_num3(k3)+1;
    end
    Rw_avg1(k1)=Rw(k1)/(W_num1(k1)+B_num1(k1)+Q_num1(k1));
    Rw_avg2(k2)=Rw(k2)/(W_num2(k2)+B_num2(k2)+Q_num2(k2));
    Rw_avg3(k3)=Rw(k3)/(W_num3(k3)+B_num3(k3)+Q_num3(k3));
    Rb_avg1(k1)=Rl(k1)/(W_num1(k1)+B_num1(k1)+Q_num1(k1));
    Rb_avg2(k2)=Rl(k2)/(W_num2(k2)+B_num2(k2)+Q_num2(k2));
    Rb_avg3(k3)=Rl(k3)/(W_num3(k3)+B_num3(k3)+Q_num3(k3));  % update the change caused by the new user 
    if Rq(i)>0   % for QoS user
        Q_chs1(i)=Rb_avg1(k1);
        B_chs1(i)=0;
        Q_chs2(i)=Rb_avg2(k2);
        B_chs2(i)=0;
        Q_chs3(i)=Rb_avg2(k3); 
        B_chs3(i)=0;
    else  % for BE user
        B_chs1(i)=Rb_avg1(k1);
        Q_chs1(i)=0;
        B_chs2(i)=Rb_avg2(k2);
        Q_chs2(i)=0;
        B_chs3(i)=Rb_avg2(k3);
        Q_chs3(i)=0;
    end
    B_min1(i)=min(Rb_avg1);
    B_min2(i)=min(Rb_avg2);
    B_min3(i)=min(Rb_avg3); 
    W_min1(i)=min(Rw_avg1);
    W_min2(i)=min(Rw_avg2);
    W_min3(i)=min(Rw_avg3);
% update the arrival of next user
    T(i)=random('exponential',lamda);  % user arriving interval follows Possion 
    T(i)=round(T(i)*10)/10; 
    w(i+1)=w(i)+T(i); i=i+1; % user arriving moment
end

%% plot the results
figure (1)
plot(B_min,'-*m')
hold on
plot(B_min1,'-or')
hold on
plot(B_min2,'-+g');
hold on
plot(B_min3,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Minimum LTE-U User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('QoS & AC','U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么

figure (2)
plot(W_min,'-*m')
hold on
plot(W_min1,'-or')
hold on
plot(W_min2,'-+g');
hold on
plot(W_min3,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Minimum WiFi User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('QoS & AC','U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么
x=get(gca,'xlim');
y=2;
plot(x,[y y],'k')  % threshold
 
figure (3)
L_chs=B_chs+Q_chs;
L_chs1=B_chs1+Q_chs1;
L_chs2=B_chs2+Q_chs2;
L_chs3=B_chs3+Q_chs3;
plot(L_chs,'-*m')
hold on
plot(L_chs1,'-or') 
hold on
plot(L_chs2,'-+g');
hold on
plot(L_chs3,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Newly Arrived LTE-U User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('QoS & AC','U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么

figure (4)
Q_chs(find(Q_chs==0))=[];
Q_chs1(find(Q_chs1==0))=[];
Q_chs2(find(Q_chs2==0))=[];
Q_chs3(find(Q_chs3==0))=[];
m=length(Q_chs);
n=length(Q_chs1);
Q_chs=[Q_chs,zeros(1,n-m)];
plot(Q_chs,'-*m')
hold on
plot(Q_chs1,'-or') 
hold on
plot(Q_chs2,'-+g');
hold on
plot(Q_chs3,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Newly Arrived LTE-U Type-1 User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('QoS & AC','U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么

figure (5)
B_chs(find(B_chs==0))=[];
B_chs1(find(B_chs1==0))=[];
B_chs2(find(B_chs2==0))=[];
B_chs3(find(B_chs3==0))=[];
m=length(B_chs);
n=length(B_chs1);
B_chs=[B_chs,zeros(1,n-m)];
plot(B_chs,'-*m')
hold on
plot(B_chs1,'-or')
hold on
plot(B_chs2,'-+g');
hold on
plot(B_chs3,'-sb');
grid on
xlabel('LTE-U user arrivals'),ylabel('Newly Arrived LTE-U Type-2 User Throughput (Mbps)')  %设置x轴和y轴的名称
legend('QoS & AC','U-CCS','RS','CS')  %区分一下蓝色和红色分别代表什么


% using utility function to measure the satisification degree of LTE-U QoS and BE users
% figure (6)
% rmax=20;
% k=1;
% C=1;



