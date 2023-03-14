% clear all;close all;clc;            %������б������ر����д��ڣ�  �������ڵ�����
% x = 0 + (100 - 0).* rand(100, 2);

function APcl = AP(x)
N=size(x,1);              %NΪ��������������������ݵ�ĸ���
M=N*N-N;                  %N�������M���������ߣ����ǵ���i��k�ʹ�k��i�ľ�������ǲ�һ����
s=zeros(M,3);             %����һ��M��3�е������,���ڴ�Ÿ������ݵ����������ƶ�

j=1;                      %ͨ��forѭ����s��ֵ����һ�б�ʾ���i���ڶ���Ϊ�յ�k��������Ϊi��k�ĸ�ŷʽ������Ϊ���ƶ�
for i=1:N
    for k=[1:i-1,i+1:N]
        s(j,1)=i;s(j,2)=k;
        %         s(j,3)=acos(1-pdist2(x(i,:),x(k,:),'cosine'));
        s(j,3)=-sum((x(i,:)-x(k,:)).^2);
        j=j+1;
    end
end
p=median(s(:,3));           %pΪ����s�����е��м�ֵ�����������ƶ�ֵ����λ��������λ����Ϊpreference,������������ʵĴصĸ���
tmp=max(max(s(:,1)),max(s(:,2)));           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=-Inf*ones(N,N);                           %-Inf������󣬶���SΪN*N�����ƶȾ��󣬳�ʼ��ÿ��ֵΪ�������
for j=1:size(s,1)                           %��forѭ����sת��ΪS��S��i��j����ʾ��i����j�����ƶ�ֵ
    S(s(j,1),s(j,2))=s(j,3);
end
nonoise=1;                                  %�˴���ѡ����������������S��i��j��=S��j,i������������ȥ���漸�д���
if length(p)==1                             %����preference
    for i=1:N
        S(i,i)=p;
    end
else
    for i=1:N
        S(i,i)=p(i);
    end
end
% Allocate space for messages ,etc
dS=diag(S);                                                 %�����������S�жԽ���Ԫ����Ϣ
A=zeros(N,N);
R=zeros(N,N);
%Execute parallel affinity propagation updates
convits=50;maxits=500;                               %���õ���������Ϊ500�Σ������������Ϊ50
e=zeros(N,convits);dn=0;i=0;                         %eѭ���ؼ�¼50�ε�����Ϣ��dn=1��Ϊһ��ѭ�������źţ�i������¼ѭ������
while ~dn
    i=i+1;
    %Compute responsibilities
    Rold=R;                                          %��Rold���¸���ǰ��R
    AS=A+S;                                          %A(i,j)+S(i,j)
    [Y,I]=max(AS,[],2);                              %���AS��ÿ�е����ֵ��ŵ�������Y�У�ÿ�����ֵ��AS�е�������ŵ�������I��
    for k=1:N
        AS(k,I(k))=-realmax;                         %��AS��ÿ�е����ֵ��Ϊ������󸡵������Ա�������Ѱ��ÿ�еĵڶ���ֵ
    end
    [Y2,I2]=max(AS,[],2);                            %���ԭAS��ÿ�еĵڶ���ֵ����Ϣ
    R=S-repmat(Y,[1,N]);                             %����R,R(i,k)=S(i,k)-max{A(i,k')+S(i,k')}��k'~=k�������������Ϊi��Ĵ����ĵ��ʺϳ̶�
    
    for k=1:N                                        %eg:��һ����AS(1,2)���,AS(1,3)�ڶ���
        R(k,I(k))=S(k,I(k))-Y2(k);                   %so R(1,1)=S(1,1)-AS(1,2); R(1,2)=S(1,2)-AS(1,3); R(1,3)=S(1,3)-AS(1,2).............
    end                                              %��������R��R��ֵ���ʾk��ô�ʺ���Ϊi �Ĵ����ģ���k�����ʺ�i�ĵ㣬��R(i,k)��ֵΪ��
    lam=0.5;
    R=(1-lam)*R+lam*Rold;                            %��������ϵ������ֹĳЩ����³��ֵ�������
    %Compute availabilities
    Aold=A;
    Rp=max(R,0);                                     %��R(k,k)�⣬��R�еĸ�����Ϊ0�����Բ��ʺ͵ĵ�Ĳ��ʺϳ̶���Ϣ
    for k=1:N
        Rp(k,k)=R(k,k);
    end
    A=repmat(sum(Rp,1),[N,1])-Rp;                    %����A(i,k),�Ƚ�ÿ�д�����Ķ�����������Ϊi~=k,����Ҫ��ȥ��ӵ�Rp(i,k)
    dA=diag(A);
    A=min(A,0);                         %��A(k,k)���⣬�����Ĵ���0��Aֵ����Ϊ0
    for k=1:N
        A(k,k)=dA(k);
    end
    A=(1-lam)*A+lam*Aold;               %��������ϵ������ֹĳЩ����³��ֵ�������
    %Check for convergence
    E=((diag(A)+diag(R))>0);
    e(:,mod(i-1,convits)+1)=E;               %��ѭ��������������E�������e�У�ע����ѭ����Ž��������һ��ѭ���ó���E�ŵ�N*50��e����ĵ�һ�У���51�εĽ���ַŵ���һ��
    K=sum(E);                                %ÿ��ֻ����������convits��ѭ��������Ա�����ж��Ƿ���������50�����Ĵؽ�������䡣%%%%%%%%%%%%%%%%
    if i>=convits || i>=maxits               %�ж�ѭ���Ƿ���ֹ
        se=sum(e,2);                         %seΪ��������E��convits�ε��������
        unconverged=(sum((se==convits)+(se==0))~=N);%���еĵ�Ҫô����50�ζ�����A+R>0��Ҫôһֱ��С���㣬��������Ϊ������
        if (~unconverged&&(K>0))||(i==maxits) %����50�β��䣬���д����Ĳ����򳬹����ѭ������ʱѭ����ֹ��
            dn=1;
        end
    end
end
I=find(diag(A+R)>0);                        %I��ÿ���ص����ģ��� %���������ѭ������ȷ��������Щ�������Ϊ�����ĵ㣬��find�����ҳ���Щ��1���ĵ�,�����demo��I=[2,4],
K=length(I); % Identify exemplars                                                                                                           %���ڶ�����͵��ĸ���Ϊ��������Ĵ�����
if K>0                                      %��������ĵĸ�������0
    [~,c]=max(S(:,I),[],2);                 %ȡ��S�еĵڶ������У����2��4�е�ÿ�е����ֵ�������һ�еڶ��е�ֵ���ڵ�һ�е����е�ֵ����˵����һ�����ǵڶ������ǹ�����
    c(I)=1:K;                               % Identify clusters c(2)=1,c(4)=2(��2����Ϊ��һ�������ģ���4����Ϊ��2��������)
    % Refine the final set of exemplars and clusters and return results
    for k=1:K
        ii=find(c==k);                      %k=1ʱ�����ֵ�1��2��3����Ϊ�����ڵ�һ����
        [y,j]=max(sum(S(ii,ii),1));         %k=1ʱ ��ȡ��S��1��2��3�к�1��2��3����ɵ�3*3�ľ��󣬷ֱ����3��֮��ȡ���ֵ��y��¼���ֵ��j��¼���ֵ���ڵ���
        I(k)=ii(j(1));                      %I=[2;4]
    end
    [tmp,c]=max(S(:,I),[],2);        %tmpΪ2��4����ÿ���������ɵ���������cΪÿ���������S������I���е�λ�ã�����ʾ���㵽�Ǹ����������
    c(I)=1:K;                        %c(2)=1;c(4)=2;
    tmpidx=I(c);                     %I=[2;4],c�е�1��2�滻��2��4�滻
    %(tmpidx-1)*N+(1:N)'                    %һ���������ֱ��ʾS(1,2),S(2,2),S(3,2),S(4,4),S(5,4),S(6,4)��S����ĵڼ���Ԫ��
    %sum(S((tmpidx-1)*N+(1:N)'))            %��S��S(1,2)+S(2,2)+S(3,2)+S(4,4)+S(5,4)+S(6,4)�ĺ�
    tmpnetsim=sum(S((tmpidx-1)*N+(1:N)'));  %�����㵽�����ĵ�һ����ʾ����ĸ�ֵ�ĺ���������ξ�����ʺ϶�
    tmpexpref=sum(dS(I));                   %dS=diag(S)��               %��ʾ���б�ѡΪ�����ĵĵ���ʺ϶�֮��
else
    tmpidx=nan*ones(N,1);  %nan Not A Number ������һ�����ݡ����ݴ���ʱ����ʵ�ʹ����о������ݵ�ȱʧ���߲���������ʱ���ǿ��Խ���Щȱʧ����Ϊnan
    tmpnetsim=nan;
    tmpexpref=nan;
end



netsim=tmpnetsim;                           %��Ӧ��ξ�����ʺ϶�
dpsim=tmpnetsim-tmpexpref;                  %���ݵ��������������ƶ�֮��
expref=tmpexpref;                           %��ʶ��ʾ����ƫ��֮��
idx=tmpidx;                                 %��¼��ÿ���������Ǹ������ĵ�������
unique(idx);

h = 1;
for i=unique(idx)'
    ii=find(idx==i);
    APcl{h} = ii';
    h = h + 1;
    
end
end

