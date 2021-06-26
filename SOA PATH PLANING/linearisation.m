%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：路径去冗余线性化函数
%%入口参数：需要去冗余线性化的路径
%%出口参数：去冗余线性化的路径,使其平直划
%%说明：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newx,newy]=linearisation(x,y)
[px,py]=size(x); %测量矩阵的行数和列数
newx=zeros(px,py);
newy=zeros(px,py);
fitness=0;
judge=1;
x3=[8 8 12 12 2 6 6 4 7 2 9 16 16 18 7 12 8 9 9 8];         %设置环境模型地图1
y3=[8 12 12 8 2 2 4 6 15 15 16 10 18 1 6 2 14 14 18 18];        %x3、y3、x4、y4 为所有障碍物边界的相邻顶点组成的线段。
x4=[8 12 12 8 6 6 2 7 2 4 16 16 9 7 12 18 9 9 8 8];
y4=[12 12 8 8 2 4 2 15 15 6 10 18 16 6 2 1 14 18 18 14];

% x3=[3 3 7 7 2 6 6 12 15 10 9 16 16 13 2 7 17 18 18 17];         %设置环境模型地图2
% y3=[8 12 12 8 14 14 16 2 11 11 18 12 20 1 6 2 2 2 6 6];        %x3、y3、x4、y4 为所有障碍物边界的相邻顶点组成的线段。
% x4=[3 7 7 3 6 6 2 15 10 12 16 16 7 2 7 13 18 18 17 17];
% y4=[12 12 8 8 14 16 14 11 11 2 12 20 18 6 2 1 2 6 6 2];

% x3=[0 3 3 0 0 8 8 0 10 13 13 10 9 16 16 9 14 15 15 14 6 10 10 6 2 6 6 2 17 20 20 17];         %设置环境模型地图三
% x4=[3 3 0 0 8 8 0 0 13 13 10 10 16 16 9 9 15 15 14 14 10 10 6 6 6 6 2 2 20 20 17 17];
% y3=[8 8 9 9 2 2 4 4 4 4 1 1 13 13 11 11 9 9 5 5 9 9 7 7 15 15 12 12 11 11 13 13];        %x3、y3、x4、y4 为所有障碍物边界的相邻顶点组成的线段。
% y4=[8 9 9 8 2 4 4 2 4 1 1 4 13 11 11 13 9 5 5 9 9 7 7 9 15 12 12 15 11 13 13 11];

for i=1:1:px
    judge=1;
    for j=1:1:py-1
        x1(j)=x(i,j);
        y1(j)=y(i,j);
    end
    for j=2:1:py
        x2(j-1)=x(i,j);
        y2(j-1)=y(i,j);
    end
%     for j=2:1:py-1
%         if sum(chack(x1(j),y1(j),x2(end),y2(end),x3,y3,x4,y4))==0   %判段该点与终点连线是否与所有障碍物边界不相交
%             for g=1:1:size(x(i,j+1:end-1),2)
%                 x(i,j+g)=rand*(x(i,end)-x(i,j))+x(i,j);
%             end
%             x(i,j+1:end-1)=sort(x(i,j+1:end-1));
%             for k=j+1:1:py-1
%                 y(i,k)=((x(i,k)-x(i,j))/(x(i,end)-x(i,j)))*(y(i,end)-y(i,j))+y(i,j);
%             end
%         end
%     end
    for j=2:1:py-1
        ch=chack(x1(judge),y1(judge),x2(j),y2(j),x3,y3,x4,y4);   %判段两点连成的线段是否与所有障碍物边界相交
        if sum(chack(x1(judge),y1(judge),x2(end),y2(end),x3,y3,x4,y4))==0   %判段该点与终点连线是否与所有障碍物边界不相交
            for g=1:1:size(x(i,judge+1:end-1),2)
                x(i,judge+g)=rand*(x(i,end)-x(i,judge))+x(i,judge);
            end
            x(i,judge+1:end-1)=sort(x(i,judge+1:end-1));
            for k=judge+1:1:py-1
                y(i,k)=((x(i,k)-x(i,judge))/(x(i,end)-x(i,judge)))*(y(i,end)-y(i,judge))+y(i,judge);
            end
        end
        fitness=sum(ch);
%         disp(fitness);
        if fitness>0
            for g=1:1:size(x(i,judge+1:j-1),2)
                x(i,judge+g)=rand*(x(i,j)-x(i,judge))+x(i,judge);
            end
            x(i,judge+1:j-1)=sort(x(i,judge+1:j-1));
            for k=judge+1:1:j-1
                y(i,k)=((x(i,k)-x(i,judge))/(x(i,j)-x(i,judge)))*(y(i,j)-y(i,judge))+y(i,judge);
            end
            judge=j;
%             disp(x(i,j));
        end
    end
    newx(i,:)=x(i,:);
    newy(i,:)=y(i,:);
end


%             x(i,judge+1:j)=sort(x(i,judge+1:j));