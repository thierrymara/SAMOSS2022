function [PCE] = SPCE(Fact,y,Degre_Init,Level_Init,Type)

[Nsample,Nvar]=size(Fact);
if nargin==2
    Degre_Init=2;
    Level_Init=2;
elseif nargin==3
    Level_Init=2;
    Nbredinconnu=factorial(Degre_Init+Level_Init)/(factorial(Degre_Init)*factorial(Level_Init));
    while Nbredinconnu<1500
        Level_Init=Level_Init+1;
        Nbredinconnu=factorial(Degre_Init+Level_Init)/(factorial(Degre_Init)*factorial(Level_Init));
    end
elseif nargin==4
    for p=1:Nvar
      Type{p}=[];
    end
end

Matdegre = PCE_Structure(Nvar,Degre_Init,Level_Init);%Create the Matrix of Polynomial Degrees

%I compute the coefficients of the PC associated with discrete variables up
%to the maximal degree
PDFHandled=['uniform','legendre','normal','gaussian','hermite','discrete'];
cnt=0;
DV=[];
for p=1:Nvar
    if findstr('discrete',lower(Type{p}))
        Fact_r=sort(Fact(:,p));
        Ind=find((Fact_r(2:end)-Fact_r(1:end-1))>0);
        DiscreteValues=Fact_r([0;Ind]+1);
        if (min(DiscreteValues)~=-1)|(max(DiscreteValues)~=1)
            warning('Discrete Variable Values should be scaled within [-1,1]!')
            PCE=[];
            return;
        end
        WeightOfValues=[];
        for n=1:length(DiscreteValues),
            Ind=find(Fact_r==DiscreteValues(n));
            WeightOfValues=[WeightOfValues;length(Ind)/length(Fact_r)];
        end
        Coeff_polynome_discrete_variable(DiscreteValues,WeightOfValues,p);
        cnt=cnt+1;
        DV(cnt).Values=DiscreteValues;
        DV(cnt).Variable=p;
        Ind1=find(Matdegre(:,p)<(length(DiscreteValues)-1));
        Matdegre=Matdegre(Ind1,:);
    elseif isempty(findstr(lower(Type{p}),PDFHandled))
        Fact(:,p)=(Fact(:,p)-mean(Fact(:,p)))/std(Fact(:,p));
    end
end

Display='On';
if strcmp(Display,'On')
    fprintf('The number of terms in the SPCE is: %d\n',size(Matdegre,1))
end
PCE.Mat=Matdegre;%Only usefull for passing the first test (while)
Matdegre=[Matdegre;Matdegre(end,:)];%Only usefull for passing the first test (while)
Degre_max=Degre_Init;%Initial max PCE degree
Level_max=Level_Init;%Initial max PCE level of interactions
Reject=[];
while size(Matdegre,1)~=size(PCE.Mat,1)%Stopping criterion
    if strcmp(Display,'On')
        fprintf('Polynomial Degree:  %d\t Interaction Level:  %d\n', Degre_max,Level_max)
    end
    Info.Step=-1;
    while size(PCE.Mat,1)~=size(Matdegre,1)
        Matdegre=PCE.Mat;
        [PCE,Rejet] = Filtre_matrix_pce(Fact,y,Matdegre,Type,Info);%
        [obs,Ind]=sort(PCE.Cr(2:end).^2,1,'descend');
        PCE.Mat=[PCE.Mat(1,:);PCE.Mat(Ind+1,:)];
        Info.Step=1;
        Reject=[Reject;Rejet];
    end
    
    Index=find(sum(PCE.Mat')>=Degre_max-1);
    if ~isempty(Index)
        Degre_max=Degre_max+2;%Increase the maximal degree of PCE
        for j=1:length(Index)
            Ind=find(PCE.Mat(Index(j),:)>0);
            for jj=1:length(Ind)
                Vect1=PCE.Mat(Index(j),:);
                Vect1(1,Ind(jj))=Vect1(1,Ind(jj))+2;
                Test=1;
                for tt=1:size(PCE.Mat,1)
                    if sum(abs(PCE.Mat(tt,:)-Vect1))==0
                        Test=-1;
                        break;
                    end
                end
                if Test==1
                    PCE.Mat=[PCE.Mat;PCE.Mat(Index(j),:)];
                    PCE.Mat(end,Ind(jj))=PCE.Mat(end,Ind(jj))+2;
                end
            end
        end
    end
    
    MatOne=PCE.Mat./max(PCE.Mat,1);
    Param=find(sum(MatOne)>0);%Relevant Inputs
    Indix=find(sum(MatOne')==Level_max);%Is the Level interactions = current max level?
    if ~isempty(Indix)
        Level_max=min(Nvar,Level_max+1);%Increase the maximal level of interactions
        for j=1:length(Indix)
            Indux=find(PCE.Mat(Indix(j),Param)==0);
            d=Degre_max-sum(PCE.Mat(Indix(j),:));
            for jj=1:d
                for jjj=1:length(Indux)
                    Vect=PCE.Mat(Indix(j),:);
                    Vect(1,Param(Indux(jjj)))=jj;
                    Test=1;
                    for tt=1:size(PCE.Mat,1)
                        if sum(abs(PCE.Mat(tt,:)-Vect))==0
                            Test=-1;
                            break;
                        end
                    end
                    if Test==1
                        PCE.Mat=[PCE.Mat;PCE.Mat(Indix(j),:)];
                        PCE.Mat(end,Param(Indux(jjj)))=jj;
                    end
                end
            end
        end
    end
    %Filters for the Discrete Variables
    for p=1:length(DV)
        Ind1=find(PCE.Mat(:,DV(p).Variable)<(length(DV(p).Values)-1));
        PCE.Mat=PCE.Mat(Ind1,:);
    end
end
Info.Step=-1;%From the 1st step (with re-ordering)
Matdegre=[PCE.Mat;Reject];
while size(PCE.Mat,1)~=size(Matdegre,1)
    Matdegre=PCE.Mat;
    [PCE] = Filtre_matrix_pce(Fact,y,Matdegre,Type,Info);
    Info.Step=1;
    Info.KIC=PCE.KIC;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SPCE,Rejet] = Filtre_matrix_pce(Fact,Sortie,Matdegre,Type,Info)

[Nsample,Nvar]=size(Fact);
Matdegre1=Matdegre;
Ordre_max=max(sum(Matdegre'));

for p=1:Nvar
    if strcmp(lower(Type{p}),'uniform')|strcmp(lower(Type{p}),'legendre');
        MatriceCoeff(p).Value = Coeff_polynome2legendre(max(max(Matdegre)));
        for d=1:Ordre_max
            Alpha1(d,p) = 1/(2*d+1);
        end
    elseif strcmp(lower(Type{p}),'gaussian')|strcmp(lower(Type{p}),'normal')|strcmp(lower(Type{p}),'hermite');
        MatriceCoeff(p).Value = Coeff_polynome2hermite(max(max(Matdegre)));
        for d=1:Ordre_max
            Alpha1(d,p) = factorial(d);
        end
    elseif findstr('discrete',lower(Type{p}))
        load(['PC_Coefficients_DV' num2str(p) '.mat'])
        MatriceCoeff(p).Value = PC_Coefficients_DV;
        for d=1:size(MatriceCoeff(p).Value,1)
            Alpha1(d,p) = 1;
        end
    else
        MatriceCoeff(p).Value = coeff_polynomeqqc_gs(Fact(:,p),max(Matdegre(:,p)));
        for d=1:Ordre_max
            Alpha1(d,p) = 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=Sortie;
if Info.Step<=0
    Variance(size(Matdegre,1),size(Sortie,2))=0;
    for r = 2:size(Matdegre,1)
        Matrice=ones(size(Sortie,1),1);
        Ind=find(Matdegre(r,:)>0);
        if ~isempty(Ind)
            for n=1:length(Ind)
                for p = 1:Matdegre(r,Ind(n))+1
                    Mat_puiss(:,p) = Fact(:,Ind(n)).^(p-1);
                end
                Vecteur = Mat_puiss*(MatriceCoeff(Ind(n)).Value(Matdegre(r,Ind(n))+1,1:Matdegre(r,Ind(n))+1)');
                Matrice = Matrice.*Vecteur;
                clear Mat_puiss
            end
        else
            disp('Erreur quelque part : voir Matdegre')
        end
        Variance(r,1)=(corr(Matrice,y).^2);
    end
    Matdegre1=Matdegre;
    
    [Variance,Indice]=sort(Variance(2:end),1,'descend');
    Indice=[1;Indice+1];
    Matdegre=Matdegre1(Indice,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y=Sortie;
    Variance(size(Matdegre,1),size(Sortie,2))=0;
    for r = 2:size(Matdegre,1)
        %waitbar(r/size(Matdegre,1),figh)
        Matrice=ones(size(Sortie,1),1);
        Ind=find(Matdegre(r,:)>0);
        if ~isempty(Ind)
            for n=1:length(Ind)
                for p = 1:Matdegre(r,Ind(n))+1
                    Mat_puiss(:,p) = Fact(:,Ind(n)).^(p-1);
                end
                Vecteur = Mat_puiss*(MatriceCoeff(Ind(n)).Value(Matdegre(r,Ind(n))+1,1:Matdegre(r,Ind(n))+1)');
                Matrice = Matrice.*Vecteur;
                clear Mat_puiss
            end
        else
            disp('Erreur quelque part : voir Matdegre')
        end
        y1=y;
        x_std=(Matrice-mean(Matrice))/std(Matrice);
        y_std=(y-mean(y))/std(y);
        Ai=polyfit(x_std,y_std,1);
        Variance(r,1)=(Ai(1).^2)*var(y)/var(Sortie);
        [obs,Coeff,y1]=Main_effect1(Matrice-mean(Matrice),y-mean(y),1);
        y=y1;
    end
    Matdegre1=Matdegre;
    
    [Variance,Indice]=sort(Variance(2:end),1,'descend');
    Variance=[0;Variance];
    Indice=[1;Indice+1];
    Matdegre1=Matdegre1(Indice,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_std=(Sortie-mean(Sortie))/std(Sortie);
Rejet=[];
Matdegre=Matdegre1(1,:);
Matrice=ones(size(Sortie,1),1);
KIC_Prec=Inf;
for r=2:size(Matdegre1,1)
    Matdegre=[Matdegre;Matdegre1(r,:)];
    Matrice=[Matrice,ones(size(Sortie,1),1)];
    Ind=find(Matdegre(end,:)>0);
    if ~isempty(Ind)
        for n=1:length(Ind)
            for p = 1:Matdegre(end,Ind(n))+1
                Mat_puiss(:,p) = Fact(:,Ind(n)).^(p-1);
            end
            Vecteur = Mat_puiss*(MatriceCoeff(Ind(n)).Value(Matdegre(end,Ind(n))+1,1:Matdegre(end,Ind(n))+1)');
            Matrice(:,end) = Matrice(:,end).*Vecteur/(sqrt(Alpha1(Matdegre(end,Ind(n)),Ind(n))));%orthonormalized
            clear Mat_puiss
        end
    else
        disp('Erreur quelque part : voir Matdegre')
    end
    Sigma2e=1;
    Sigma2e_prec=0;
    counter=0;
    while max(abs(Sigma2e_prec-Sigma2e)./Sigma2e)>0.0001
        Sigma2e_prec=Sigma2e;
        Cr1=((Matrice')*Matrice)\((Matrice')*y_std);%OLS Estimate
        Sigma2e = sum((y_std-Matrice*Cr1).^2)/Nsample;
        counter=counter+1;
        if counter==100
            break;
        end
    end
    
    DetJ=Inf;
    N=0;
    while or(isinf(DetJ),DetJ==0)
        DetJ=det((Matrice')*Matrice/(Nsample^N*Sigma2e));%I divide by Nsample^N to avoid DetJ=Inf. This is corrected below in the KIC
        if isinf(DetJ)
            N=N+1;
        elseif DetJ==0
            N=N-1/20;
        end
    end
    KIC=Nsample*log(Sigma2e)+N*size(Matrice,2)*log(Nsample)-size(Matrice,2)*log(2*pi)+log(abs(DetJ));%Kashyap IC, as small as possible

    if isinf(KIC),KIC=-Inf;Result=Cr1.^2,Sigma2e,abs(DetJ),disp('là');break;end
    
    if (KIC<KIC_Prec)&(Cr1(end)^2>Cr1(1)^2/10)%%%%%%%%%%%%%%%%%%%%%%
        KIC_Prec=KIC;
    else
        Matrice=Matrice(:,1:end-1);
        Rejet=[Rejet;Matdegre(end,:)];
        Matdegre=Matdegre(1:end-1,:);
        Cr1=((Matrice')*Matrice)\((Matrice')*y_std);%OLS Estimate
    end
end
Sigma2e=mean((y_std-Matrice*Cr1).^2);
Result=Cr1.^2;
Result(1) = 0;%Remove the mean value
SPCE(1).Cr=Cr1;
SPCE(1).Mat=Matdegre;
SPCE(1).Var=sum(Result)*var(Sortie);
SPCE(1).Res=sum((y_std-Matrice*Cr1).^2)/(Nsample*sum(Result));
SPCE(1).Muy=mean(Sortie);
SPCE(1).Sigmay=std(Sortie);
SPCE(1).Cov=inv((Matrice')*Matrice/Sigma2e);
SPCE(1).Sigma2e=Sigma2e;
SPCE(1).Ndraw=Nsample;
%Test of Normality
SPCE(1).Residual=(y_std-Matrice*Cr1);
SPCE(1).KIC=KIC_Prec;

