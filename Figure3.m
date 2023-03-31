
% set_path_to_MUA_data

%%%% MUA CTOA %%%%

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
    CTOA = CTOA(tmp)
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    [mu1, sig1, betas1,beta_fast1,beta_main1, graphhandles1, x_1,y_1] = reciprobit2(RT/1000);
    [res_x, idx_of_result] = knee_pt(x_1,y_1,0);
    tmp1 = find(RT>-1000/x_1(idx_of_result));
    CTOA = CTOA(tmp1)
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
end

Matr_CondT_CTOA=reshape(Matr_CondT,kkk*48,num_bins,801);

figure;
%axes('NextPlot', 'add');
color1 = [0, 48, 143]./255;      % Set colors according to your definition of
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    color = color1 * m(k) + color2 * (1 - m(k));
    plot(squeeze(mean(Matr_CondT_CTOA(:,k,:),1)), 'Color', color,'LineWidth',2); hold all
end

A=squeeze(mean(Matr_CondT_CTOA(:,:,200:400),3));
for i=1:10, 
    
    error(i)=std(A(:,i))/sqrt(size(A,1));
end


figure;
%axes('NextPlot', 'add');
color1 = [0, 48, 143]./255;      % Set colors according to your definition of
color2 = [142, 223, 255]./255;    % "light blue" and "dark blue"
m      = linspace(0, 1, 10);
ctoa = 1000:(3200 -1000)/10:3200 - (3200 -1000)/10
for k = 1:10
    color = color1 * m(k) + color2 * (1 - m(k));
   plot(ctoa(k),mean(A(:,k),1),'square','Color',color,'MarkerSize',20,'MarkerFaceColor',color); hold all,errorbar(ctoa(k),mean(A(:,k),1),error(i),'Color',color,'LineWidth',2);;hold all
end

figure, histogram(CTOA_tot,20)

%%%% MUA RT %%%%

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
end

Matr_CondT_RT=reshape(Matr_CondT,kkk*48,num_bins,801);

figure;

color1 = [0, 48, 143]./255;       
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    color = color1 * m(k) + color2 * (1 - m(k));
    plot(squeeze(mean(Matr_CondT_RT(:,k,:),1)), 'Color', color,'LineWidth',2); hold all
end

A=squeeze(mean(Matr_CondT_RT(:,:,200:400),3));
for i=1:10, 
    
    error(i)=std(A(:,i))/sqrt(size(A,1));
end


figure;
color1 = [0, 48, 143]./255;      
color2 = [142, 223, 255]./255;    
m      = linspace(0, 1, 10);
rt = min(RT_tot):(max(RT_tot)-min(RT_tot))/10:max(RT_tot) - (max(RT_tot)-min(RT_tot))/10
for k = 1:10
    color = color1 * m(k) + color2 * (1 - m(k));
   plot(rt(k),mean(A(:,k),1),'square','Color',color,'MarkerSize',20,'MarkerFaceColor',color); hold all,errorbar(rt(k),mean(A(:,k),1),error(i),'Color',color,'LineWidth',2);;hold all
end

figure, histogram(RT_tot,10)


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
end

Matr_CondT_TA=reshape(Matr_CondT,kkk*48,num_bins,801);

figure;

color1 = [0, 48, 143]./255;       
color2 = [142, 223, 255]./255;
m      = linspace(0, 1, num_bins);
for k = 1:num_bins
    color = color1 * m(k) + color2 * (1 - m(k));
    plot(squeeze(mean(Matr_CondT_TA(:,k,:),1)), 'Color', color,'LineWidth',2); hold all
end

A=squeeze(mean(Matr_CondT_TA(:,:,200:400),3));
for i=1:10, 
    
    error(i)=std(A(:,i))/sqrt(size(A,1));
end


figure;
color1 = [0, 48, 143]./255;      
color2 = [142, 223, 255]./255;    
m      = linspace(0, 1, 10);
ta = min(TA_tot):(max(TA_tot)-min(TA_tot))/10:max(TA_tot) - (max(TA_tot)-min(TA_tot))/10
for k = 1:10
    color = color1 * m(k) + color2 * (1 - m(k));
   plot(ta(k),mean(A(:,k),1),'square','Color',color,'MarkerSize',20,'MarkerFaceColor',color); hold all,errorbar(ta(k),mean(A(:,k),1),error(i),'Color',color,'LineWidth',2);;hold all
end

figure, histogram(TA_tot,10)


