
function [cg_curve,zbestx,zbesty]=SOA(N,Max_iteration,lb,ub,dim,Px0,Py0,Px1,Py1)
format long
sizepop=N;%��Ⱥ��ģ
maxgen=Max_iteration;%����������
m=dim;%�ռ�ά��
Umax=0.9500;%���������ֵ
Umin=0.0111;%��С������ֵ
Wmax=0.9;%���Ȩ��
Wmin=0.1;%��СȨ��
popmax=ub;
popmin=lb;
%��ʼ����Ⱥ����
[x,y]=Path_init(sizepop,m,Px0,Py0,Px1,Py1);   %������ʼ��Ⱥ
[x,y]=linearisation(x,y);%·��ȥ����
%[x,y]=[a,b];
fitness=calfitvalue(x,y);%��Ӧ�ȼ���

[bestfitness,bestindex]=max(fitness);%Ѱ�Ҿ��������Ӧ�ȵĸ���  ����Сֵʱ[bestfitness bestindex]=min(fitness);
zbestx=x(bestindex,:);%ȫ�����
zbesty=y(bestindex,:);%ȫ�����
gbestx=x;%�������
gbesty=y;%�������
fitnesszbest=bestfitness;%ȫ�������Ӧ��ֵ
fitnessgbest=fitness;%���������Ӧ��ֵ
%��һ��Ϊ��ʼ��������ʽ����
%����Ѱ��
Dix=zeros(sizepop,m);%0*rand(sizepop,m)
Diy=zeros(sizepop,m);%0*rand(sizepop,m)
Dix(1,:)=1;%��һ��ȫΪ1
Diy(1,:)=1;%��һ��ȫΪ1
step_lengthx=zeros(sizepop,m);%����������ʼ��
step_lengthy=zeros(sizepop,m);%����������ʼ��
Diegox=zeros(sizepop,m);%��������
Diegoy=zeros(sizepop,m);%��������
Dialtx=zeros(sizepop,m);%�����������
Dialty=zeros(sizepop,m);%�����������
Diprox=zeros(sizepop,m);%Ԥ������
Diproy=zeros(sizepop,m);%Ԥ������
cg_curve=zeros(maxgen,1);
cg_curve(1)=fitnesszbest;

for t=2:maxgen%��ʽ������ʼ��Ϊ��һ���ǳ�ʼ�����Ե�����2��ʼ
     [x,y]=linearisation(x,y);%·��ȥ����
     for i=1:sizepop%����ÿһ������
         W=Wmax-t*(Wmax-Wmin)/maxgen;%����Ȩֵ
         Diegox(i,2:(dim-1))=sign(gbestx(i,2:(dim-1))-x(i,2:(dim-1)));%�����������
         Diegoy(i,2:(dim-1))=sign(gbesty(i,2:(dim-1))-y(i,2:(dim-1)));%�����������
         Dialtx(i,2:(dim-1))=sign(zbestx(2:(dim-1))-x(i,2:(dim-1)));%�����������
         Dialty(i,2:(dim-1))=sign(zbesty(2:(dim-1))-y(i,2:(dim-1)));%�����������
         if calfitvalue(gbestx(i,:),gbesty(i,:))<=calfitvalue(x(i,:),y(i,:))%Ԥ���������,����Сֵ����>=,���ֵ<=
            Diprox(i,2:(dim-1))=-Dix(i,2:(dim-1));
            Diproy(i,2:(dim-1))=-Diy(i,2:(dim-1));
         else
            Diprox(i,2:(dim-1))=Dix(i,2:(dim-1));
            Diproy(i,2:(dim-1))=Diy(i,2:(dim-1));
         end
         Dix(i,2:(dim-1))=sign(W*Diprox(i,2:(dim-1))+rand*Diegox(i,2:(dim-1))+rand*Dialtx(i,2:(dim-1)));%�����������
         Diy(i,2:(dim-1))=sign(W*Diproy(i,2:(dim-1))+rand*Diegoy(i,2:(dim-1))+rand*Dialty(i,2:(dim-1)));%�����������
         %^ȷ�������ݶȷ���
        [Orderfitnessgbest,indexfitnessgbest]=sort(fitnessgbest,'descend');%�������и�����Ӧ�ȣ���������,�����������͸�ֵ��Ӧ��ԭλ��
        
        u=Umax-(sizepop-indexfitnessgbest(i))*(Umax-Umin)/(sizepop-i);%u��ȷ�����������������ȡ����������Ϊ0.95����СΪ0.011���õ�һ��λ�õ�������֮���������[u,1]��������u
        U=u+(1-u)*rand;%�������[u��1]��������u
        
        H(t)=(maxgen-t)/maxgen;%����������Ȩ�صı仯
        Cx(i,:)=H(t)*abs(zbestx-5*rands(1,m));%ȷ����˹�����Ĳ���
        Cy(i,:)=H(t)*abs(zbesty-5*rands(1,m));%ȷ����˹�����Ĳ���
        T=sqrt(-log(U));%%%������ȡ������-ln��u��
        
        step_lengthx(i,2:(dim-1))=Cx(i,2:(dim-1))*T;%ȷ�����������Ĵ�С
        step_lengthy(i,2:(dim-1))=Cy(i,2:(dim-1))*T;%ȷ�����������Ĵ�С
        step_lengthx(i,find(step_lengthx(i,2:(dim-1))>3*max(Cx(i,2:(dim-1)))))=3*max(Cx(i,2:(dim-1)));%���������д���������˹�����Ⱥ������ֵ�ò���ȫȡ������˹�����Ⱥ������ֵ%%�޲���
        step_lengthx(i,find(step_lengthx(i,2:(dim-1))<0))=0;%����������С��0�Ĳ���ȫȡ0%%�޲���
        step_lengthy(i,find(step_lengthy(i,2:(dim-1))>3*max(Cy(i,2:(dim-1)))))=3*max(Cy(i,2:(dim-1)));%���������д���������˹�����Ⱥ������ֵ�ò���ȫȡ������˹�����Ⱥ������ֵ%%�޲���
        step_lengthy(i,find(step_lengthy(i,2:(dim-1))<0))=0;%����������С��0�Ĳ���ȫȡ0%%�޲���
        
        %����λ��
        x(i,2:(dim-1))=x(i,2:(dim-1))+Dix(i,2:(dim-1)).*step_lengthx(i,2:(dim-1));%�����˸���
        y(i,2:(dim-1))=y(i,2:(dim-1))+Diy(i,2:(dim-1)).*step_lengthy(i,2:(dim-1));%�����˸���
        x(i,find(x(i,2:(dim-1))>popmax))=popmax;%���ڱ�������ȡ����
        x(i,find(x(i,2:(dim-1))<popmin))=popmin;%С�ڱ�������ȡ����
        y(i,find(y(i,2:(dim-1))>popmax))=popmax;%���ڱ�������ȡ����
        y(i,find(y(i,2:(dim-1))<popmin))=popmin;%С�ڱ�������ȡ����
        
        %[x(i,:),y(i,:)]=linearisation(x(i,:),y(i,:));
        fitness(i)=calfitvalue(x(i,:),y(i,:));%������Ӧ��ֵ
        %�������Ÿ���
        if fitness(i)>fitnessgbest(i) && fitness(i)>0
            gbestx(i,:)=x(i,:); %%%%%%%%%%%%%%%
            gbesty(i,:)=y(i,:);
            fitnessgbest(i) = fitness(i);
        end
        %Ⱥ�����Ÿ���
        if fitness(i) > fitnesszbest && fitness(i)>0
            zbestx = x(i,:);
            zbesty = y(i,:);
            fitnesszbest = fitness(i);
        end
    cg_curve(t)=fitnesszbest;
         
    end
    
end

end