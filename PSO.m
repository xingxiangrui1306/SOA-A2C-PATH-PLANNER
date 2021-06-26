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

Vmax=1;%����ٶ�
Vmay=1;%����ٶ�
noP=N;%������
wMax=0.9;%������
wMin=0.2;%��С����
c1=2;%ѧϰ����
c2=2;%ѧϰ����

% Initializations
iter=Max_iteration;%����������
%vel=zeros(noP,dim);%��ʼ���������ڴ�����Ӻͱ���
velx=zeros(noP,dim);%��ʼ���������ڴ�����Ӻ�x����
vely=zeros(noP,dim);%��ʼ���������ڴ�����Ӻ�y����
pBestScore=zeros(noP,1);%����nop��nop�е�0����
pBestx=zeros(noP,dim);%��ʼ�������Ÿ������ӵ��������λ��
pBesty=zeros(noP,dim);%��ʼ�������Ÿ������ӵ��������λ��
gBestx=zeros(1,dim);%��ʼ�������Ÿ���Ⱥ���λ��
gBesty=zeros(1,dim);%��ʼ�������Ÿ���Ⱥ���λ��
cg_curve=zeros(1,iter);%�������ڷ���ÿһ�ε��������Ž�
bestx=zeros(iter,dim);
besty=zeros(iter,dim);

% Random initialization for agents.��ʼ��·��
%pos=initialization(noP,dim,ub,lb); 
[a,b]=Path_init(noP,dim,Px0,Py0,Px1,Py1);   %������ʼ��Ⱥ
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
     x(i,:)=(x(i,:).*(~(Flagxub+Flagxlb)))+ub.*Flagxub+lb.*Flagxlb;%�����������޵�ֵ����ֵ����������           
     Flagyub=y(i,:)>ub;
     Flagylb=y(i,:)<lb;
     y(i,:)=(y(i,:).*(~(Flagyub+Flagylb)))+ub.*Flagyub+lb.*Flagylb;%�����������޵�ֵ����ֵ����������
         
        %Calculate objective function for each particle��������ӵ���Ӧ��
        fitness=calfitvalue(x(i,:),y(i,:));

        if(pBestScore(i)<fitness && fitness>0)%���¸����ӵ�����ֵ
            pBestScore(i)=fitness;
            pBestx(i,:)=x(i,:);
            pBesty(i,:)=y(i,:);
        end
        if(gBestScore<fitness && fitness>0)%������Ⱥ����
            gBestScore=fitness;
            gBestx=x(i,:);
            gBesty=y(i,:);
        end
    end
%     display(gBestScore);
%     h=plot(gBestx,gBesty,'r-');hold on
%     legend(h) 
    %Update the W of PSO
    w=wMax-l*((wMax-wMin)/iter);%���¶�̬����Ȩֵ�����Եݼ�
    %Update the Velocity and Position of particles�������ӵ��ٶȺ�λ��
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