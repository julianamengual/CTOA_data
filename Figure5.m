%%%% CTOA %%%%
d=dir([cd '/' '*data.mat']);
kkk=0;
num_bins = 10
for iii=[1:length(d)]
    kkk=kkk+1;
    load(d(iii).name,'data*','Behavior');
    iii
    beh=Behavior(:,1);
    CTOA = Behavior(:,3);
    RT = Behavior(:,4);
    tmp=find(beh==1 & CTOA<3000 & CTOA > 1000);
    RT=RT(tmp);
    CTOA = CTOA(tmp);
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    [mu1, sig1, betas1,beta_fast1,beta_main1, graphhandles1, x_1,y_1] = reciprobit2(RT/1000);
    [res_x, idx_of_result] = knee_pt(x_1,y_1,0);
    tmp1 = find(RT>-1000/x_1(idx_of_result));
    CTOA = CTOA(tmp1);
    if iii ==1
        CTOA_tot = CTOA;
    else
        CTOA_tot = cat(1,CTOA_tot,CTOA);
    end
    dataT=dataT(:,:,tmp1);
    dataC=dataC(:,:,tmp1);
    for bin = 1:num_bins
        trial{bin}=find(CTOA>prctile(CTOA,ceil((100/num_bins)*(bin-1))) & CTOA<prctile(CTOA,ceil((100/num_bins)*(bin))));
        for chan=1:48
            Matr_CondT(kkk,chan,bin,:)=squeeze(mean(dataT(chan,:,trial{bin}),3));
        end
    end
    
    for chan=1:48
        for i=1:num_bins,
            x=trial{i};
            x=x(randperm(length(x)));
            A(chan,i,:,1:length(x))=dataT(chan,:,x);
        end
    end
    firingRatesT{kkk}=A;
    clear A B %
end

for k=1:length(firingRatesT);
    k
    A=firingRatesT{k};
    firingRates_CTOA(k,:,:,:,1:size(A,4))=A;    
end

Matr_CondT_CTOA=reshape(Matr_CondT,kkk*48,num_bins,801);
s=size(firingRates_CTOA);
firingRates_CTOA=reshape(firingRates_CTOA,kkk*48,s(3),s(4),s(5));
combinedParams = {{1, [1 2]}, {2}};
margNames = {'CTOA',  'C-Indep'};
margColours = [23 100 171; 187 20 25; ]/256;

[W,V,whichMarg] = dpca(Matr_CondT_CTOA, 100, ...
    'combinedParams', combinedParams);
explVar = dpca_explainedVariance(Matr_CondT_CTOA, W, V, ...
    'combinedParams', combinedParams);
optimalLambda = dpca_optimizeLambda(Matr_CondT_CTOA(:,:,1:10:end), firingRates_CTOA(:,:,1:10:end,:), size(firingRates_CTOA,4), ...
    'combinedParams', combinedParams, ...
    'simultaneous', true, ...
    'numRep', 10, ... 
    'filename', []);
% close all
Cnoise = dpca_getNoiseCovariance(Matr_CondT_CTOA(:,:,1:10:end), ...
    firingRates_CTOA(:,:,1:10:end,:), size(firingRates_CTOA,4), 'simultaneous', true);

 S=10
 decodingClasses = {[(1:S)']};
[accuracy_CTOA,~,TPR] = dpca_classificationAccuracy(Matr_CondT_CTOA(:,:,1:10:end), firingRates_CTOA(:,:,1:10:end,:), size(firingRates_CTOA,4), ...
    'numComps',3,...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', true, ...
    'numRep', 100, ...        % increase to 100
    'filename', []);

%% Plot
[a b] = sort(explVar.margVar(1,:),'descend')
W_ctoa = W(:,b(1:3));

X = Matr_CondT_CTOA(:,:);
X = bsxfun(@minus, X, mean(X,2));
k=0;
    for j=[1:3]
        k=k+1;
       P(k,:)=X'*W_ctoa(:,j);
    end
    
 PP=reshape(P,3,num_bins,801);
 

figure;

color1 = [0, 48, 143]./255;     
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    subplot(2,3,[1:3])
    color = color1 * m(k) + color2 * (1 - m(k));
    plot3(squeeze(PP(1,k,:)),squeeze(PP(2,k,:)),squeeze(PP(3,k,:)), 'Color', color,'LineWidth',2.5); hold all
    plot3(squeeze(PP(1,k,400)),squeeze(PP(2,k,400)),squeeze(PP(3,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    plot3(squeeze(PP(1,k,800)),squeeze(PP(2,k,800)),squeeze(PP(3,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    xlabel('dPCA1')
    ylabel('dPCA2')
    zlabel('dPCA3')
    grid on;
    
    for i=1:3
        subplot(2,3,i+3);
        plot(-400:400,squeeze(PP(i,k,:)), 'Color', color,'LineWidth',2.5);hold all; ylabel('dPCA1')
        plot(0,squeeze(PP(i,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;ylabel(['dPCA ' num2str(i)])
        plot(400,squeeze(PP(i,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;
    end
end


%%%% RT %%%%


d=dir([cd '/' '*data.mat']);
kkk=0;
num_bins = 10
for iii=[1:length(d)]
    kkk=kkk+1;
    load(d(iii).name,'data*','Behavior');
    iii
    beh=Behavior(:,1);
    RT = Behavior(:,4);
    tmp=find(beh==1);
    RT=RT(tmp);
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    [mu1, sig1, betas1,beta_fast1,beta_main1, graphhandles1, x_1,y_1] = reciprobit2(RT/1000);
    [res_x, idx_of_result] = knee_pt(x_1,y_1,0);
    tmp1 = find(RT>-1000/x_1(idx_of_result));
    RT = RT(tmp1);
    if iii ==1
        RT_tot = RT;
    else
        RT_tot = cat(1,RT_tot,RT);
    end
    dataT=dataT(:,:,tmp1);
    dataC=dataC(:,:,tmp1);
    for bin = 1:num_bins
        trial{bin}=find(RT>prctile(RT,ceil((100/num_bins)*(bin-1))) & RT<prctile(RT,ceil((100/num_bins)*(bin))));
        for chan=1:48
            Matr_CondT(kkk,chan,bin,:)=squeeze(mean(dataT(chan,:,trial{bin}),3));
        end
    end
         for chan=1:48
        for i=1:num_bins,
            x=trial{i};
            x=x(randperm(length(x)));
            A(chan,i,:,1:length(x))=dataT(chan,:,x);
        end
    end
    firingRatesT{kkk}=A;
    clear A B %
end


Matr_CondT_RT=reshape(Matr_CondT,kkk*48,num_bins,801);

for k=1:length(firingRatesT);
    k
    A=firingRatesT{k};
    firingRates_RT(k,:,:,:,1:size(A,4))=A;    
end

Matr_CondT_RT=reshape(Matr_CondT,kkk*48,num_bins,801);
s=size(firingRates_RT);
firingRates_RT=reshape(firingRates_RT,kkk*48,s(3),s(4),s(5));
combinedParams = {{1, [1 2]}, {2}};
margNames = {'RT',  'C-Indep'};
margColours = [23 100 171; 187 20 25; ]/256;

[W,V,whichMarg] = dpca(Matr_CondT_RT, 100, ...
    'combinedParams', combinedParams);
explVar = dpca_explainedVariance(Matr_CondT_RT, W, V, ...
    'combinedParams', combinedParams);
optimalLambda = dpca_optimizeLambda(Matr_CondT_RT(:,:,1:10:end), firingRates_RT(:,:,1:10:end,:), size(firingRates_RT,4), ...
    'combinedParams', combinedParams, ...
    'simultaneous', true, ...
    'numRep', 10, ... 
    'filename', []);
% close all
Cnoise = dpca_getNoiseCovariance(Matr_CondT_RT(:,:,1:10:end), ...
    firingRates_RT(:,:,1:10:end,:), size(firingRates_RT,4), 'simultaneous', true);

 S=10
 decodingClasses = {[(1:S)']};
[accuracy_RT,~,TPR] = dpca_classificationAccuracy(Matr_CondT_RT(:,:,1:10:end), firingRates_RT(:,:,1:10:end,:), size(firingRates_RT,4), ...
    'numComps',3,...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', true, ...
    'numRep', 100, ...        % increase to 100
    'filename', []);

%% Plot
[a b] = sort(explVar.margVar(1,:),'descend')
W_rt = W(:,b(1:3));

X = Matr_CondT_RT(:,:);
X = bsxfun(@minus, X, mean(X,2));
k=0;
    for j=[1:3]
        k=k+1;
       P(k,:)=X'*W_rt(:,j);
    end
    
 PP=reshape(P,3,num_bins,801);
 

figure;

color1 = [0, 48, 143]./255;     
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    subplot(2,3,[1:3])
    color = color1 * m(k) + color2 * (1 - m(k));
    plot3(squeeze(PP(1,k,:)),squeeze(PP(2,k,:)),squeeze(PP(3,k,:)), 'Color', color,'LineWidth',2.5); hold all
    plot3(squeeze(PP(1,k,400)),squeeze(PP(2,k,400)),squeeze(PP(3,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    plot3(squeeze(PP(1,k,800)),squeeze(PP(2,k,800)),squeeze(PP(3,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    xlabel('dPCA1')
    ylabel('dPCA2')
    zlabel('dPCA3')
    grid on;
    
    for i=1:3
        subplot(2,3,i+3);
        plot(-400:400,squeeze(PP(i,k,:)), 'Color', color,'LineWidth',2.5);hold all; ylabel('dPCA1')
        plot(0,squeeze(PP(i,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;ylabel(['dPCA ' num2str(i)])
        plot(400,squeeze(PP(i,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;
    end
end

%% TA %%

%%%% MUA TA %%%%

d=dir([cd '/' '*data.mat']);
kkk=0;
num_bins = 10
for iii=[1:length(d)]
    kkk=kkk+1;
    load(d(iii).name,'data*','Behavior');
    iii
    beh=Behavior(:,1);
    TA = Behavior(:,9);
    RT = Behavior(:,4);
    TA = TA.*14*sqrt(2)/2; % TA to degrees
    tmp=find(beh==1);
    TA=TA(tmp);
    RT = RT(tmp);
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    [mu1, sig1, betas1,beta_fast1,beta_main1, graphhandles1, x_1,y_1] = reciprobit2(RT/1000);
    [res_x, idx_of_result] = knee_pt(x_1,y_1,0);
    tmp1 = find(RT>-1000/x_1(idx_of_result) & TA <20);
    TA = TA(tmp1);
    if iii ==1
        TA_tot = TA;
    else
        TA_tot = cat(1,TA_tot,TA);
    end
    dataT=dataT(:,:,tmp1);
    dataC=dataC(:,:,tmp1);
    for bin = 1:num_bins
        trial{bin}=find(TA>prctile(TA,ceil((100/num_bins)*(bin-1))) & TA<prctile(TA,ceil((100/num_bins)*(bin))));
        for chan=1:48
            Matr_CondT(kkk,chan,bin,:)=squeeze(mean(dataT(chan,:,trial{bin}),3));
        end
    end
              for chan=1:48
        for i=1:num_bins,
            x=trial{i};
            x=x(randperm(length(x)));
            A(chan,i,:,1:length(x))=dataT(chan,:,x);
        end
    end
    firingRatesT{kkk}=A;
    clear A B %
end

Matr_CondT_TA=reshape(Matr_CondT,kkk*48,num_bins,801);

for k=1:length(firingRatesT);
    k
    A=firingRatesT{k};
    firingRates_TA(k,:,:,:,1:size(A,4))=A;    
end

Matr_CondT_TA=reshape(Matr_CondT,kkk*48,num_bins,801);
s=size(firingRates_TA);
firingRates_TA=reshape(firingRates_TA,kkk*48,s(3),s(4),s(5));
combinedParams = {{1, [1 2]}, {2}};
margNames = {'TA',  'C-Indep'};
margColours = [23 100 171; 187 20 25; ]/256;

[W,V,whichMarg] = dpca(Matr_CondT_TA, 100, ...
    'combinedParams', combinedParams);
explVar = dpca_explainedVariance(Matr_CondT_TA, W, V, ...
    'combinedParams', combinedParams);
optimalLambda = dpca_optimizeLambda(Matr_CondT_TA(:,:,1:10:end), firingRates_TA(:,:,1:10:end,:), size(firingRates_TA,4), ...
    'combinedParams', combinedParams, ...
    'simultaneous', true, ...
    'numRep', 10, ... 
    'filename', []);
% close all
Cnoise = dpca_getNoiseCovariance(Matr_CondT_TA(:,:,1:10:end), ...
    firingRates_TA(:,:,1:10:end,:), size(firingRates_TA,4), 'simultaneous', true);

 S=10
 decodingClasses = {[(1:S)']};
[accuracy_TA,~,TPR] = dpca_classificationAccuracy(Matr_CondT_TA(:,:,1:10:end), firingRates_TA(:,:,1:10:end,:), size(firingRates_TA,4), ...
    'numComps',3,...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', true, ...
    'numRep', 100, ...        % increase to 100
    'filename', []);

%% Plot
[a b] = sort(explVar.margVar(1,:),'descend')
W_rt = W(:,b(1:3));

X = Matr_CondT_TA(:,:);
X = bsxfun(@minus, X, mean(X,2));
k=0;
    for j=[1:3]
        k=k+1;
       P(k,:)=X'*W_rt(:,j);
    end
    
 PP=reshape(P,3,num_bins,801);
 

figure;

color1 = [0, 48, 143]./255;     
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    subplot(2,3,[1:3])
    color = color1 * m(k) + color2 * (1 - m(k));
    plot3(squeeze(PP(1,k,:)),squeeze(PP(2,k,:)),squeeze(PP(3,k,:)), 'Color', color,'LineWidth',2.5); hold all
    plot3(squeeze(PP(1,k,400)),squeeze(PP(2,k,400)),squeeze(PP(3,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    plot3(squeeze(PP(1,k,800)),squeeze(PP(2,k,800)),squeeze(PP(3,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all
    xlabel('dPCA1')
    ylabel('dPCA2')
    zlabel('dPCA3')
    grid on;
    
    for i=1:3
        subplot(2,3,i+3);
        plot(-400:400,squeeze(PP(i,k,:)), 'Color', color,'LineWidth',2.5);hold all; ylabel('dPCA1')
        plot(0,squeeze(PP(i,k,400)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;ylabel(['dPCA ' num2str(i)])
        plot(400,squeeze(PP(i,k,800)),'square', 'Color', color,'Markersize',10,'MarkerFaceColor',color); hold all;
    end
end



