clc,clear;
disp('Length scale of Station spacing is 1');

%% configure network parameters
num=4;  % Number of Stations
n=num+1;  % Number of Stations and an AP
x1=0;y1=0;  % AP location
r=3;  % Range of Radio
motion_scale=100;  % Time Scale of random motion
frame_size=8;  % Average Transmission Time (packet size)
max_simutime=1000;  % Simulation Time in ms (SlotTime=1)
disp('Performing calculations -- this may take a while, so please be patient.')
ph=0:0.02:2*pi;
figure
hold on
title(' DCF ');
j=sqrt(-1);
for i=1:n
    h2(i)=plot(exp(j*ph)*(1+r)); % blue node
    h4(i)=plot(exp(j*ph)*(1+r)); % handle for the range ellipse (cyan) - when node is transmitting
    h5(i)=plot(exp(j*ph)*(1+r),'LineWidth',3,'color','r'); % Data packet - RED colour
    h6(i)=plot(exp(j*ph)*(1+r),'LineWidth',3,'color','g'); % ACK packet  - GREEN 
end

%% definitions and initialization 
alpha=0.5^(1/motion_scale);
Eb2=0.001;
beta=(1-alpha^2-Eb2)/(2*Eb2-Eb2*alpha+(1-alpha^2)*alpha);
ER2=(1-beta^2)*Eb2;
traffic=1.0;  % traffic=1,saturated
ACK_length = 4;
SIFS=4;  % SIFS
SlotTime=1;  % slot time(us)
DIFS=SIFS+4*SlotTime;  % DIFS
PCWmin=15;  % minimum contention window
PCWmax=1023;  % maximum contention window
PCW=(log(PCWmin+1)/log (2))*ones(1,n);  %backoff stage
range=r*ones(1,n);  
state=zeros(1,n);  % initial state of all stations
prev_state=5*ones(1,n);  
current_frame_length=zeros(1,n);  % Sifs+Difs+rand(1,n)*frame_size; 
Timer_DIFS=zeros(1,n);
Timer_SIFS=zeros(1,n);
Timer_ACK=zeros(1,n);  % timer
BC=zeros(1,n);  % backoff counter 
BC_Init=zeros(1,n);  % backoff counter 
pos_x=randn(1,num)+x1;  
pos_y=randn(1,num)+y1;  
pos_x_change=sqrt(Eb2)*randn(1,num);
pos_y_change=sqrt(Eb2)*rand(1,num); % position & position change of stations
pos_x=[pos_x,x1];
pos_y=[pos_y,y1];
pos_x_change=[pos_x_change,0];
pos_y_change=[pos_y_change,0];  % position & position change of stations and AP
within=ones(n,n);  % whether a pair of stations could communicate
total_transmissions = 0;
successful_transmission = 0;
total_collisions = 0;
unreachable_packets = 0;
total_acks = 0;
ack_collisions = 0;
unreachable_acks = 0;
successful_acks = 0;  
total_packets=0;  % record transmission rsults

%% run
for counter=1:max_simutime  
    t0=clock;  % current time    
    pos_x_change(1,1:num)=beta*pos_x_change(1,1:num)+randn(1,num)*sqrt(ER2);
    pos_y_change(1,1:num)=beta*pos_y_change(1,1:num)+randn(1,num)*sqrt(ER2);
%     pos_x_change=[pos_x_change,0];
%     pos_y_change=[pos_y_change,0];
    pos_x=pos_x*alpha+pos_x_change;
    pos_y=pos_y*alpha+pos_y_change;  % update locations of STAs
    for i = 1:n
        for k = 1:n        
            if((pos_x(i)-pos_x(k))^2 + (pos_y(i)-pos_y(k))^2 > range(i)^2)                     
                within(i,k) = 0;            
            end
        end
    end  % update within information
   
%% State - 0 idle (nothing to send)   
    temp=rand(1,n);
    transit0to1=zeros(1,n);
    transit0to1(state==0 & temp<traffic)=1;
    temp=floor(rand(1,n).*(2.^PCW-1));  % select back off counter from 0 to minimum contention window
    BC(transit0to1>0 & state==0) = temp(transit0to1>0 & state==0)*SlotTime;
    BC_Init(transit0to1>0 & state==0) = temp(transit0to1>0 & state==0)*SlotTime;
    
%% State - 1 Medium Sensing     
   sending=zeros(1,n);
   sending(state>=4 | state<=-1 )=1;
   Busy_media=sending*within;  % Busy_media determines which nodes are within range and sending   
   transit1to2 = zeros(1,n) ;    
   transit1to2(state == 1 & Busy_media < 1) = 1; % when medium idle, transit to state 2  
   
%% State - 2 Difs Timer
   Timer_DIFS(transit1to2>0) = DIFS; % Wait for "DIFS" amount of time before sending the Data packet. 
   transit2to1 = zeros(1,n);
   transit2to1(state==2 & Busy_media>0)=1; % If sensing busy, back to state 1      
   Timer_DIFS(state==2) = Timer_DIFS(state==2)-1;  % Counting down for DIFS amount of time
   % print the DIFS counter
   if (state==2)
       if (Timer_DIFS == DIFS-1)
           disp('wait DIFS before initial transmit');
       end
       disp(Timer_DIFS(state==2));
   end
   transit2to3 = zeros(1,n);
   transit2to3(state== 2 & Timer_DIFS==0) = 1; % when DIFS=0, transit to state 3
     
 %% State - 3 Backoff Counter
%    BC(state==3 & Busy_media<1) = BC(state==3 & Busy_media<1)-1; % back off
   % print the backoff counter
   if (state==3 & Busy_media<1)
      if BC(state==3 & Busy_media<1)==BC_Init(state==3 & Busy_media<1)
        disp('ALL stations Decrement their backoff counter');
      end
      disp(BC(state==3 & Busy_media<1));
   end
   transit3to1=zeros(1,n);
   transit3to1(state==3 & Busy_media>0)=1; % If sensing busy, back to state 1       
   transit3to5 = zeros(1,n);
   transit3to5(state==3 & BC==0) = 1; % when BC<0, transit to state 5
   BC(state==3 & Busy_media<1) = BC(state==3 & Busy_media<1)-1; % back off
   if any(transit3to5>0)    
        temp1=rand(1,n)*frame_size*2;
        current_frame_length=temp1*1000/130; % frame length
        initial_cfl=current_frame_length; % initial frame length, no change in run
        initial_cfl1=temp1;
%         temp1=ceil(rand(1,n)*n); % select a destination randomly
        current_frame_dest = n*ones(1,n);  % the des must be AP
        for i=1:n
            while (current_frame_dest(i)==i | current_frame_dest(i)==0)
                current_frame_dest(i)=ceil(rand(1,1)*n); % avoid no des or itself being des
           end
        end
   end
   
  %% State - 5 transmission
   for i=1:n  
       for j=1:n
           if (state(i)==5 & state(j)==5 & j~=i) % collision occurs
               state(j)=6;state(i)=6; % collision state
               if prev_state(i)==5
                   disp([' OOPS COLLISION OCCURS: station ' num2str(j) ' is also transmitting']);
               else
                   disp([' OOPS COLLISION OCCURS: station ' num2str(i) ' is also transmitting']);
               end
           end
       end
       if (state(i)==5 & within(i,current_frame_dest(i))<1)
           state(i)=7; % the destination node is out of range
           disp(['station ' num2str(current_frame_dest(i)) ' is out of ranges from station ' num2str(i)]);
       end
   end
   for i=1:n
       if(state(i) == 5 & current_frame_length(i) > 0)
           if current_frame_length(i)==initial_cfl(i)
               disp(['station ' num2str(i) ' sending data to station ' num2str(current_frame_dest(i))]) % tx  
           end
           disp(current_frame_length(i));
       end
   end
   current_frame_length(state==5 |state==6 | state==7) = current_frame_length(state==5 | state==6 | state==7)-1; % one slot passed
   transit5to0 = zeros(1,n);
   transit5to0(state ==5 & current_frame_length < 0) = 1; % sending successfully
   transit5to4 = zeros(1,n);
   transit6to0 = zeros(1,n);
   transit7to0 = zeros(1,n);
   PCW(state == 6) = PCW(state == 6)+1;
   transit6to0(state == 6) = 1;
   transit7to0(state == 7) = 1;
   for i = 1:n
       if(transit5to0(i) > 0)
%            if(state(i)==6)
%                PCW(state == 6) = PCW(state == 6)+1; % collision, backoff stage increase by 1
%                transit6to0((i)) = 1;
%            elseif(state(i)==7)
%                transit7to0((i)) = 1;
%            else
               total_packets=total_packets+initial_cfl1(i);
               transit5to4(current_frame_dest(i)) = 1; % des is going to send ack
               Timer_SIFS(current_frame_dest(i)) = SIFS+1;  
               Timer_ACK(current_frame_dest(i)) = ACK_length;
               ACK_dest(current_frame_dest(i)) = i; % wait for ack
%            end
       end
   end
   for i=1:n
       if (state(i)==-2 | state(i)==6)
           cont_size=2.^PCW(i)-1; % collision, disp increased contention window
           disp(['contention window of station ' num2str(i) ' increases to ' num2str(cont_size) ])
       end 
   end
   maxpcw=(log(PCWmax+1)/log (2))*ones(1,n);
   if PCW(state == -2 | state ==6)>= maxpcw(state == -2 | state ==6)
       PCW(state == -2 | state ==6)= maxpcw(state == -2 | state ==6);
       disp ('maximum contention window reached'); % reach the maximum contention window
   end
   total_transmissions = total_transmissions + length(state(state >= 5 & transit5to0 > 0)); % tx number
   successful_transmission = successful_transmission + length(state(state == 5 & transit5to0 >0)); % succ tx number
   total_collisions = total_collisions + length(state(state == 6 & transit5to0 > 0)); % collision number
   unreachable_packets = unreachable_packets + length(state(state == 7 & transit5to0 > 0)); % unreachable number
   state(transit5to4 > 0) = 4; % if tx suc, des trans to state 4, sending ack
   Timer_SIFS(state == 4 & transit5to4 > 0) = SIFS;     
   Timer_SIFS(state == 4) = Timer_SIFS(state == 4) - 1;  
   % print the SIFS count down   
   for j=1:n
       if (state(j)==4 & Timer_SIFS(j)>=-1)
           if (Timer_SIFS(j)==SIFS-1)
               disp(['station ' num2str(j) ' wait SIFS before sending ACK']);
           end
           disp(Timer_SIFS(state==4 & Timer_SIFS(j)>=0));
       end
   end
   transit4tom1 = zeros(1,n);
   transit4tom1(state == 4 & Timer_SIFS <= 0) = 1; % when SIFS finished, send ack, trans to state -1
     
  %% State - -1 sending ack
   for i = 1:n           
       if (state(i) == -1 && Busy_media(ACK_dest(i)) > 1)
           state(i) = -2;
       end
       if (state(i) == -1 && within(i,ACK_dest(i)) < 1)
           state(i) = -3;
       end
   end    
   Timer_ACK(state <= -1) = Timer_ACK(state <= -1) - 1;
   for j=1:n
       if (state(j)==-1 & Timer_SIFS(j)<=0)
           if Timer_ACK(state <= -1) ==ACK_length-1
               disp(['station ' num2str(j) ' sending ACK to ' num2str(ACK_dest(j)) ]);
           end
           disp(Timer_ACK(state==-1 & Timer_SIFS<=0));
       end
   end
   transitm1to0 = zeros(1,n);
   transitm1to0(state <= -1 & Timer_ACK <=0) = 1;
   total_acks = total_acks + length(state(state <= -1 & transitm1to0 > 0));
   successful_acks = successful_acks + length(state(state == -1 & transitm1to0 > 0));
   ack_collisions = ack_collisions + length(state(state == -2 & transitm1to0 > 0));
   unreachable_acks = unreachable_acks + length(state(state == -3 & transitm1to0 > 0));
%    total_packets=total_packets+initial_cfl(state==5 & transit5to0>0);
  
 %% States Conversion   
   state(transit0to1 > 0)  =  1;
   state(transit3to1 > 0)  =  1;
   state(transit1to2 > 0)  =  2;
   state(transit2to3 > 0)  =  3;
   state(transit3to5 > 0)  =  5;     
   state(transit2to1 > 0)  =  1;
   state(transit5to0 > 0)  =  0;
   state(transit6to0 > 0)  =  0;
   state(transit7to0 > 0)  =  0;
   state(transit5to4 > 0)  =  4;
   state(transit4tom1 > 0) = -1;
   state(transitm1to0 > 0) =  0;    
   prev_state=state;
%    state
   SHOW_DCF
   drawnow  
   while etime(clock,t0)<0.1
   end     
end

total_transmissions
successful_transmission
total_collisions
unreachable_packets
total_acks
successful_acks
ack_collisions
unreachable_acks 
total_packets
throughput=total_packets/max_simutime*1000  % 130Mbit/s channel data rate
avg_throughput=throughput/num
   %STATE TABLE 
%-3: unreachable ACK
%-2: ACK collision
%-1: sending an ACK
%0: idle (nothing to send)
%1: Waiting for free media (have something to send)
%2: Free media detected, wait for DIFS
%3: Random back-off
%4: Waiting for SIFS to send ACK
%5: sending successfully
%6: sending to a busy station
%7: sending to an out-of-range node
  