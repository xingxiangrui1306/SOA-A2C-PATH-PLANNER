 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：初始化族群
%%入口参数：种群数量  基因数量
%%出口参数：初始种群
%%说明：
    %%初始种群的个点的X轴坐标与Y轴坐标分开存放，分别放在矩阵 x,y中，作为函数返回值返回
    %%初始种群的产生，除去起始点与终止点两点，其他点的x轴、y轴随机产生，并从大到小进行排列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=Path_init(popsize,chromlength,x0,y0,xn,yn)
i=1;
for j=1:1:popsize %种群数量从一开始初始化，每个个体就是一条路径
    while i<=chromlength %基因数量从1开始，每一个基因代表一个路径坐标点
        if i==1
            x(j,i)=x0; %每个种群的第一个基因都是起点
            y(j,i)=y0;
            i=i+1;
        elseif i<chromlength-1 %中间的n-2个基因
            if rand<0.5
                x(j,i)=(xn-x(j,i-1))*rand+x(j,i-1);%x不断靠近目标
                y(j,i)=(yn-y0)*rand; %y随机搜索整个起点到终点
            else
                y(j,i)=(yn-y(j,i-1))*rand+y(j,i-1);
                x(j,i)=(xn-x0)*rand; 
            end
            fit=calfitvalue([x(j,i-1),x(j,i)],[y(j,i-1),y(j,i)]); %进行一次适应度检测
            if fit==0
%                 i=i-1;
            else
                i=i+1;
            end
        elseif i==chromlength-1
            if rand<0.5
                x(j,i)=(xn-x(j,i-1))*rand+x(j,i-1);
                y(j,i)=(yn-y0)*rand;
            else
                y(j,i)=(yn-y(j,i-1))*rand+y(j,i-1);
                x(j,i)=(xn-x0)*rand;
            end
            fit=calfitvalue([x(j,i-1),x(j,i),xn],[y(j,i-1),y(j,i),yn]); %进行一次适应度检测
            if fit==0
                 i=1;%这里导致了所有初始化一定不经过图形
            else
                i=i+1;
            end
        elseif i==chromlength %最后一个
            x(j,i)=xn;
            y(j,i)=yn;
            i=i+1;
        end
    end
    i=1;
end
% x=20.0*rand(popsize,chromlength);
% y=20.0*rand(popsize,chromlength);
% x(:,1)=0;
% y(:,1)=0;
% x(:,chromlength)=20;
% y(:,chromlength)=20;
% [px,py]=size(x);
% % if rand<0.5
% %     for i=1:1:px
% %     x(i,:)=sort(x(i,:));
% % %     y(i,:)=sort(y(i,:));
% %     end 
% % else
% %     for i=1:1:px
% % %     x(i,:)=sort(x(i,:));
% %     y(i,:)=sort(y(i,:));
% %     end 
% % end
% for i=1:1:px
% x(i,:)=sort(x(i,:));
% y(i,:)=sort(y(i,:));
% end