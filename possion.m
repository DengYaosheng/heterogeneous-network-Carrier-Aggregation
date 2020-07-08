clc,clear,clf
fidin=fopen('D:\Program Files (x86)\MATLAB\MATLAB Compiler Runtime\v714\bin\win32\Possion_data.txt');
fidout=fopen('Possion.txt','w');
while (~feof(fidin))
    tline=fgetl(fidin);
    if(numel(tline)==0)
       continue
    end
    if(double(tline(1))>=48 && double(tline(1))<=57)
        fprintf(fidout,'%s\n\n',tline);
        continue
    end
end
fclose(fidout);
test=importdata('Possion.txt');%存放读入的possion数据
% count=zeros(1,100);
% for k=1:100
% for j=1:33
%     for i=1:100
%         if test(i,j)<=0.1*k && test(i,j)>(k-1)*0.1
%             count(k)=count(k)+1;
%         end
%     end
% end
% end
% for i=2:100
%     count(i)=count(i)+count(i-1);
% end
% for i=1:100
%     
% E(i)=count(i)./100./(0.1*i);
% end

count1=zeros(length(test(:,1)),100);%用来存放每隔0.1时间的计数值
for i=1:length(test(:,1))
for k=1:100
for j=1:length(test(1,:))
        if (test(i,j)<=0.1*k) && (test(i,j)>0)
           count1(i,k)=count1(i,k)+1;  
    end
end
end
end

for i=1:100
    m(i)=mean(count1(:,i));%同一时刻的均值
    v(i)=var(count1(:,i));%同一时刻的方差
end

i=0.1:0.1:10;
plot(i,m,'gx');
axis([0,10,0,21]);
hold on
y=polyfit(i,m,1)%最小二乘法线性拟合求lamda
y1=polyval(y,i);
plot(i,y1,'k')

plot(i,v,'.r');
y=polyfit(i,v,1)%最小二乘法线性拟合求lamda
y2=polyval(y,i);
plot(i,y2,'b')

Tn=zeros(length(test(:,1)),length(test(1,:)));%用来存放前后两次计数的时间间隔
Tn(:,1)=test(:,1);
for i=2:length(test(1,:))
    Tn(:,i)=test(:,i)-test(:,i-1);
end
N=0;S=0;
for i=1:length(test(:,1))
   for j=1:length(test(1,:))
       if Tn(i,j)>0
           S=S+Tn(i,j);
           N=N+1;
       end
   end
end
lamda=S\N%用指数分布求lamda
    

            