 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������ƣ���ʼ����Ⱥ
%%��ڲ�������Ⱥ����  ��������
%%���ڲ�������ʼ��Ⱥ
%%˵����
    %%��ʼ��Ⱥ�ĸ����X��������Y������ֿ���ţ��ֱ���ھ��� x,y�У���Ϊ��������ֵ����
    %%��ʼ��Ⱥ�Ĳ�������ȥ��ʼ������ֹ�����㣬�������x�ᡢy��������������Ӵ�С��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=Path_init(popsize,chromlength,x0,y0,xn,yn)
i=1;
for j=1:1:popsize %��Ⱥ������һ��ʼ��ʼ����ÿ���������һ��·��
    while i<=chromlength %����������1��ʼ��ÿһ���������һ��·�������
        if i==1
            x(j,i)=x0; %ÿ����Ⱥ�ĵ�һ�����������
            y(j,i)=y0;
            i=i+1;
        elseif i<chromlength-1 %�м��n-2������
            if rand<0.5
                x(j,i)=(xn-x(j,i-1))*rand+x(j,i-1);%x���Ͽ���Ŀ��
                y(j,i)=(yn-y0)*rand; %y�������������㵽�յ�
            else
                y(j,i)=(yn-y(j,i-1))*rand+y(j,i-1);
                x(j,i)=(xn-x0)*rand; 
            end
            fit=calfitvalue([x(j,i-1),x(j,i)],[y(j,i-1),y(j,i)]); %����һ����Ӧ�ȼ��
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
            fit=calfitvalue([x(j,i-1),x(j,i),xn],[y(j,i-1),y(j,i),yn]); %����һ����Ӧ�ȼ��
            if fit==0
                 i=1;%���ﵼ�������г�ʼ��һ��������ͼ��
            else
                i=i+1;
            end
        elseif i==chromlength %���һ��
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