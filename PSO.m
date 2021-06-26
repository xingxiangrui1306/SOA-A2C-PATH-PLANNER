%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Particle Swarm Optimization
function [cg_curve,bestx,besty]=PSO(N,Max_iteration,lb,ub,dim,Px0,Py0,Px1,Py1)

%PSO Infotmation

Vmax=1;%最大速度
Vmay=1;%最大速度
noP=N;%粒子数
wMax=0.9;%最大惯量
wMin=0.2;%最小惯量
c1=2;%学习因子
c2=2;%学习因子

% Initializations
iter=Max_iteration;%最大迭代次数
%vel=zeros(noP,dim);%初始化矩阵用于存放粒子和变量
velx=zeros(noP,dim);%初始化矩阵用于存放粒子和x变量
vely=zeros(noP,dim);%初始化矩阵用于存放粒子和y变量
pBestScore=zeros(noP,1);%生成nop行nop列的0矩阵
pBestx=zeros(noP,dim);%初始化矩阵存放个个粒子的自身最好位置
pBesty=zeros(noP,dim);%初始化矩阵存放个个粒子的自身最好位置
gBestx=zeros(1,dim);%初始化矩阵存放个种群最好位置
gBesty=zeros(1,dim);%初始化矩阵存放个种群最好位置
cg_curve=zeros(1,iter);%后面用于返回每一次迭代的最优解
bestx=zeros(iter,dim);
besty=zeros(iter,dim);

% Random initialization for agents.初始化路径
%pos=initialization(noP,dim,ub,lb); 
[a,b]=Path_init(noP,dim,Px0,Py0,Px1,Py1);   %产生初始种群
[x,y]=linearisation(a,b);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:noP
    pBestScore(i)=-inf;
end

% Initialize gBestScore for a minimization problem
 gBestScore=-inf;
     
    
for l=1:iter 
   
    for i=1:size(x,1)
    % Return back the particles that go beyond the boundaries of the search
    % space
     Flagxub=x(i,:)>ub;
     Flagxlb=x(i,:)<lb;
     x(i,:)=(x(i,:).*(~(Flagxub+Flagxlb)))+ub.*Flagxub+lb.*Flagxlb;%将超出上下限的值变量值换成上下限           
     Flagyub=y(i,:)>ub;
     Flagylb=y(i,:)<lb;
     y(i,:)=(y(i,:).*(~(Flagyub+Flagylb)))+ub.*Flagyub+lb.*Flagylb;%将超出上下限的值变量值换成上下限
         
        %Calculate objective function for each particle计算该粒子的适应度
        fitness=calfitvalue(x(i,:),y(i,:));

        if(pBestScore(i)<fitness && fitness>0)%跟新该粒子的最优值
            pBestScore(i)=fitness;
            pBestx(i,:)=x(i,:);
            pBesty(i,:)=y(i,:);
        end
        if(gBestScore<fitness && fitness>0)%跟新种群最优
            gBestScore=fitness;
            gBestx=x(i,:);
            gBesty=y(i,:);
        end
    end
%     display(gBestScore);
%     h=plot(gBestx,gBesty,'r-');hold on
%     legend(h) 
    %Update the W of PSO
    w=wMax-l*((wMax-wMin)/iter);%跟新动态惯性权值，线性递减
    %Update the Velocity and Position of particles更新粒子的速度和位置
    for i=1:size(x,1)
        for j=2:1:size(x,2)-1       
            velx(i,j)=w*velx(i,j)+c1*rand()*(pBestx(i,j)-x(i,j))+c2*rand()*(gBestx(j)-x(i,j));
            vely(i,j)=w*vely(i,j)+c1*rand()*(pBesty(i,j)-x(i,j))+c2*rand()*(gBesty(j)-x(i,j));
            
            if(velx(i,j)>Vmax)
                velx(i,j)=Vmax;
            end
            if(velx(i,j)<-Vmax)
                velx(i,j)=-Vmax;
            end
            if(vely(i,j)>Vmay)
                vely(i,j)=Vmay;
            end
            if(vely(i,j)<-Vmay)
                vely(i,j)=-Vmay;
            end   
            x(i,j)=x(i,j)+velx(i,j);
            y(i,j)=y(i,j)+vely(i,j);
        end
    end
    cg_curve(l)=gBestScore;
    bestx(l,:)=gBestx;
    besty(l,:)=gBesty;
%     display(gBestx);
end

end