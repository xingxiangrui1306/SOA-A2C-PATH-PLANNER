
function [cg_curve,zbestx,zbesty]=SOA(N,Max_iteration,lb,ub,dim,Px0,Py0,Px1,Py1)
format long
sizepop=N;%种群规模
maxgen=Max_iteration;%最大迭代次数
m=dim;%空间维数
Umax=0.9500;%最大隶属度值
Umin=0.0111;%最小隶属度值
Wmax=0.9;%最大权重
Wmin=0.1;%最小权重
popmax=ub;
popmin=lb;
%初始化种群个体
[x,y]=Path_init(sizepop,m,Px0,Py0,Px1,Py1);   %产生初始种群
[x,y]=linearisation(x,y);%路径去冗余
%[x,y]=[a,b];
fitness=calfitvalue(x,y);%适应度计算

[bestfitness,bestindex]=max(fitness);%寻找具有最好适应度的个体  求最小值时[bestfitness bestindex]=min(fitness);
zbestx=x(bestindex,:);%全局最佳
zbesty=y(bestindex,:);%全局最佳
gbestx=x;%个体最佳
gbesty=y;%个体最佳
fitnesszbest=bestfitness;%全局最佳适应度值
fitnessgbest=fitness;%个体最佳适应度值
%第一次为初始化，启发式搜索
%迭代寻优
Dix=zeros(sizepop,m);%0*rand(sizepop,m)
Diy=zeros(sizepop,m);%0*rand(sizepop,m)
Dix(1,:)=1;%第一行全为1
Diy(1,:)=1;%第一行全为1
step_lengthx=zeros(sizepop,m);%搜索步长初始化
step_lengthy=zeros(sizepop,m);%搜索步长初始化
Diegox=zeros(sizepop,m);%利己方向
Diegoy=zeros(sizepop,m);%利己方向
Dialtx=zeros(sizepop,m);%利他方向更新
Dialty=zeros(sizepop,m);%利他方向更新
Diprox=zeros(sizepop,m);%预动方向
Diproy=zeros(sizepop,m);%预动方向
cg_curve=zeros(maxgen,1);
cg_curve(1)=fitnesszbest;

for t=2:maxgen%正式迭代开始因为第一次是初始化所以迭代从2开始
     [x,y]=linearisation(x,y);%路径去冗余
     for i=1:sizepop%遍历每一个个体
         W=Wmax-t*(Wmax-Wmin)/maxgen;%惯性权值
         Diegox(i,2:(dim-1))=sign(gbestx(i,2:(dim-1))-x(i,2:(dim-1)));%利己方向更新
         Diegoy(i,2:(dim-1))=sign(gbesty(i,2:(dim-1))-y(i,2:(dim-1)));%利己方向更新
         Dialtx(i,2:(dim-1))=sign(zbestx(2:(dim-1))-x(i,2:(dim-1)));%利他方向更新
         Dialty(i,2:(dim-1))=sign(zbesty(2:(dim-1))-y(i,2:(dim-1)));%利他方向更新
         if calfitvalue(gbestx(i,:),gbesty(i,:))<=calfitvalue(x(i,:),y(i,:))%预动方向更新,找最小值是用>=,最大值<=
            Diprox(i,2:(dim-1))=-Dix(i,2:(dim-1));
            Diproy(i,2:(dim-1))=-Diy(i,2:(dim-1));
         else
            Diprox(i,2:(dim-1))=Dix(i,2:(dim-1));
            Diproy(i,2:(dim-1))=Diy(i,2:(dim-1));
         end
         Dix(i,2:(dim-1))=sign(W*Diprox(i,2:(dim-1))+rand*Diegox(i,2:(dim-1))+rand*Dialtx(i,2:(dim-1)));%搜索方向更新
         Diy(i,2:(dim-1))=sign(W*Diproy(i,2:(dim-1))+rand*Diegoy(i,2:(dim-1))+rand*Dialty(i,2:(dim-1)));%搜索方向更新
         %^确定经验梯度方向
        [Orderfitnessgbest,indexfitnessgbest]=sort(fitnessgbest,'descend');%排列所有个体适应度，降序排列,返回排序结果和该值对应的原位置
        
        u=Umax-(sizepop-indexfitnessgbest(i))*(Umax-Umin)/(sizepop-i);%u的确定可利用线性隶属度。最大隶属度为0.95，最小为0.011。得到一个位置的隶属度之后，随机生成[u,1]的数代替u
        U=u+(1-u)*rand;%随机生成[u，1]的数代替u
        
        H(t)=(maxgen-t)/maxgen;%迭代过程中权重的变化
        Cx(i,:)=H(t)*abs(zbestx-5*rands(1,m));%确定高斯函数的参数
        Cy(i,:)=H(t)*abs(zbesty-5*rands(1,m));%确定高斯函数的参数
        T=sqrt(-log(U));%%%这里求取根号下-ln（u）
        
        step_lengthx(i,2:(dim-1))=Cx(i,2:(dim-1))*T;%确定搜索步长的大小
        step_lengthy(i,2:(dim-1))=Cy(i,2:(dim-1))*T;%确定搜索步长的大小
        step_lengthx(i,find(step_lengthx(i,2:(dim-1))>3*max(Cx(i,2:(dim-1)))))=3*max(Cx(i,2:(dim-1)));%步长中所有大于三倍高斯隶属度函数最大值得步长全取三倍高斯隶属度函数最大值%%限步长
        step_lengthx(i,find(step_lengthx(i,2:(dim-1))<0))=0;%步长中所有小于0的步长全取0%%限步长
        step_lengthy(i,find(step_lengthy(i,2:(dim-1))>3*max(Cy(i,2:(dim-1)))))=3*max(Cy(i,2:(dim-1)));%步长中所有大于三倍高斯隶属度函数最大值得步长全取三倍高斯隶属度函数最大值%%限步长
        step_lengthy(i,find(step_lengthy(i,2:(dim-1))<0))=0;%步长中所有小于0的步长全取0%%限步长
        
        %更新位置
        x(i,2:(dim-1))=x(i,2:(dim-1))+Dix(i,2:(dim-1)).*step_lengthx(i,2:(dim-1));%个体人更新
        y(i,2:(dim-1))=y(i,2:(dim-1))+Diy(i,2:(dim-1)).*step_lengthy(i,2:(dim-1));%个体人更新
        x(i,find(x(i,2:(dim-1))>popmax))=popmax;%大于变量上限取上限
        x(i,find(x(i,2:(dim-1))<popmin))=popmin;%小于变量下限取下限
        y(i,find(y(i,2:(dim-1))>popmax))=popmax;%大于变量上限取上限
        y(i,find(y(i,2:(dim-1))<popmin))=popmin;%小于变量下限取下限
        
        %[x(i,:),y(i,:)]=linearisation(x(i,:),y(i,:));
        fitness(i)=calfitvalue(x(i,:),y(i,:));%计算适应度值
        %个体最优更新
        if fitness(i)>fitnessgbest(i) && fitness(i)>0
            gbestx(i,:)=x(i,:); %%%%%%%%%%%%%%%
            gbesty(i,:)=y(i,:);
            fitnessgbest(i) = fitness(i);
        end
        %群体最优更新
        if fitness(i) > fitnesszbest && fitness(i)>0
            zbestx = x(i,:);
            zbesty = y(i,:);
            fitnesszbest = fitness(i);
        end
    cg_curve(t)=fitnesszbest;
         
    end
    
end

end