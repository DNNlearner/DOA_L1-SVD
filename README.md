# DOA_L1-SVD
LI-SVD arithmetic in the application of DOA
M=8;%阵元数
K=3;%信源数
theta=[25,60,-45];%信号来向
snap=100;%快拍数
lamda=2;%波长
d=lamda/2;%阵元间距
fc=[2*10^3,3*10^3,2.8*10^3];%信号频率
fs=20*10^3;%采样频率
SNR=0;%信噪比
S=10.^(SNR/20)*[randn(K,snap)+1j*2*pi*fc'*(1:snap)/fs];%信号S
for i=1:K
    for n=1:M
        A(n,i)=exp(-1j*2*pi*(n-1)*d/lamda*sin(theta(i)/180*pi));
    end
end
X=A*S;
Y=awgn(X,SNR,'measured');

%% l1-SVD算法数据预处理
search_doa=-90:90;
len=length(search_doa);%设定方向网格181个方向
for j=1:len
    A_dict(:,j)=exp(-1j*2*pi*(0:(M-1))*d/lamda*sin(search_doa(j)/180*pi));%生成字典方向矩阵
end
Dk1=eye(K);
Dk2=zeros(snap-K,K);
Dk=[Dk1;Dk2];
[U,sigma,V]=svd(Y);
Ysv=Y*V*Dk;
%未知向量为Ssv1

%% 采用matlab凸优化工具箱CVX
sumvec=ones(len,1);
cvx_begin quiet
variables p q %常数
variable r(len) %矢量
variable Ssv1(len,K) complex
expressions Xsv(len,1) 
expressions Zk(M,K) complex
minimize(p+2*q)
subject to
    Zk=Ysv-A_dict*Ssv1;
    Zvec=vec(Zk);%返回包含变量所有元素的列向量
    %norm(Zvec)<=p;
    %sum(Zvec.^2)<=p;
    test=norm(Zvec);
    pow_pos(test,2)<=p;
    sumvec'*r<=q;
    for i=1:len
      Xsv(i)=norm(Ssv1(i,:));
    end
    for i=1:len
      Xsv(i)<=r(i);
    end
cvx_end
%power1=10*log10(abs(Ssv1(:,1))/max(abs(Ssv1(:,1))));%扩展0-1之间的数值从0-（-inf）可以明显看出谱峰
power1=Xsv;
plot(search_doa,power1,'r');
grid minor;
xlabel('DOA/degree');
ylabel('Power/dB');
