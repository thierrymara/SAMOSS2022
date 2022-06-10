# Octave 7.1.0, Mon Jun 06 13:59:10 2022 GMT <unknown@LAPTOP-JIV56EFT>
# cleaned up session log from JRC Summer School 2022
pkg load signal
pkg load statistics
# Moment Independent
n=1024;
% constant conditional mean
x=sobolpoints(n,3);
m=5;
y=x(:,1).*(x(:,2).^m-(1-m/(m+1)))+.1*x(:,3);
subplot(2,2,1)
cosi(x(:,1),y,4,'Constant conditional mean')
y=x(:,1)-2*x(:,1).*(x(:,2)>.6)+sqrt(x(:,3));
subplot(2,2,2)
cosi(x(:,1),y,4,'Bifurcation')
subplot(2,2,4)
y=((x(:,2)+x(:,3)).*betainv(x(:,1),x(:,2),x(:,3))-x(:,2))./sqrt(x(:,2).*x(:,3)./...
(x(:,2)+x(:,3)+1));
cosi(x(:,2),y,4,'Conditional Skewness')

[bi,di,ki,ei]=betaKS3(x,y,32,'Gfx')
deltamim(x,y,16,'skew')

cosi(x,y,4,'Conditional Skewness')

n=1024;
% constant conditional mean
x=sobolpoints(n,3);
m=5;
y=x(:,1).*(x(:,2).^m-(1-m/(m+1)))+.1*x(:,3);
subplot(2,2,1)
cosi(x(:,1),y,4,'Constant conditional mean')

deltamim(x,y,32,'density')
deltamim(x,y,struct('PartitionSize',32,'KSLevel',0),'density')

clear all
ishigami
k=3;
x=trafo(rand(2000,k));y=model(x);
clf
copdist(x,y);

x=trafo(sobolpoints(256,k));y=model(x);
copdist(x,y);
x=trafo(sobolpoints(128,k));y=model(x);
copdist(x,y);

copulasi(x,y,[100 32],'gfx')
colormap(jet)
x=trafo(sobolpoints(512,k));y=model(x);
copulasi(x,y,[100 32],'xx')
# Octave 7.1.0, Wed Jun 08 12:59:00 2022 GMT <unknown@LAPTOP-JIV56EFT>
# screening
model99
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
pause
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
hold on

[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')
[m,s,ms]=morrisEE(k,10,model,trafo,.1);
plot(abs(m),s,'*',ms,s,'o')

x=rand(1000,99); # use SA as screening
y=model(x);
Si=cosi(x,y,8);
clf
bar(Si)
Si=cosi(x(:,1:10),y,8) # just analying an input subset

# Thierry discussing HSIC
for i=1:99,Hi(i)=hsic(x(:,i),y); end
for i=1:99,Hi(i,:)=hsic(x(:,i),y); end
bar(Hi)

