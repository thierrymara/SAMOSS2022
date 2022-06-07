if exist ('OCTAVE_VERSION', 'builtin')
    warning('off')
    pkg load statistics
    %pkg load gnuplot
    struct_levels_to_print(0)
end
close all
clear all
Here=pwd;
addpath([Here,'\Benchmarks'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=8;%Number of Model Inputs;
for p=1:d
    Type{p}='Uniform';%'Gaussian';%
end
PDF.Coeff=[0.10,0.016,0.0;7.71,1,0.0;63070,100000,115600;6.95,0.0167,0.0;63,116,0.0;6.6,0.033,0.0;1120,1680,0.0;9855,12045,0.0];
PDF.Type={'Normal';'LogNormal';'Triangular';'LogNormal';'Uniform';'LogNormal';'Uniform';'Uniform'};

Nsample=1000;%Sample size
U = rand(Nsample,d);

%Transformation
X(:,1)=sqrt(2)*erfinv(2*U(:,1)-1)*PDF.Coeff(1,2)+PDF.Coeff(1,1);%Normal
%
X(:,2)=exp(sqrt(2)*erfinv(2*U(:,2)-1)*PDF.Coeff(2,2)+PDF.Coeff(2,1));%LogNormal
%
Ind1=find(U(:,3)<(PDF.Coeff(3,2)-PDF.Coeff(3,1))/(PDF.Coeff(3,3)-PDF.Coeff(3,1)));Ind2=find(U(:,3)>=(PDF.Coeff(3,2)-PDF.Coeff(3,1))/(PDF.Coeff(3,3)-PDF.Coeff(3,1)));
X(:,3)=zeros(Nsample,1);
X(Ind1,3)=PDF.Coeff(3,1)+sqrt(U(Ind1,3)*(PDF.Coeff(3,2)-PDF.Coeff(3,1))*(PDF.Coeff(3,3)-PDF.Coeff(3,1)));
X(Ind2,3)=PDF.Coeff(3,3)-sqrt((PDF.Coeff(3,3)-PDF.Coeff(3,2))*(1-U(Ind2,3))*(PDF.Coeff(3,3)-PDF.Coeff(3,1)));%Triangle
%
X(:,4)=exp(sqrt(2)*erfinv(2*U(:,4)-1)*PDF.Coeff(4,2)+PDF.Coeff(4,1));%LogNormal
%
X(:,5)=U(:,5)*(PDF.Coeff(5,2)-PDF.Coeff(5,1))+PDF.Coeff(5,1);%Uniform
%
X(:,6)=exp(sqrt(2)*erfinv(2*U(:,6)-1)*PDF.Coeff(6,2)+PDF.Coeff(6,1));%LogNormal
%
X(:,7)=U(:,7)*(PDF.Coeff(7,2)-PDF.Coeff(7,1))+PDF.Coeff(7,1);%Uniform
%
X(:,8)=U(:,8)*(PDF.Coeff(8,2)-PDF.Coeff(8,1))+PDF.Coeff(8,1);%Uniform

for k=1:Nsample
    y(k,1)=borehole(X(k,:));
end

figure;hist(y)
fprintf('The mean value of the model response is:  %5.4f\n', mean(y))
fprintf('The variance of the model response is approximately:  %5.4f\n', var(y))
