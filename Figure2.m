path_data = cd('CTOA2/Data/MUA/')
file = dir('*data.mat');

%% Each file contain two matrixes:
% Behavior: Contains behavior information of each trial (hit & misses)
 
% dataC: Smooth firing rates locked to the Cue (sensors x time x trials),
% dataT: Smooth firing rates locked to the Targets (sensors x time x trials)

%%% FIGURE 2B, MODULATION INDEX
k_chan = 0
sel_chan = 0
for session=[1:length(file)]
    load(file(session).name);
    beh = Behavior(:,1);  % 1:hit, 2:miss,
    Position = Behavior(:,2);
    % Selecting hit trial per each position
    a1 = find(beh==1 &  Position==1);
    a2 = find(beh==1 &  Position==2);
    a3 = find(beh==1 &  Position==3);
    a4 = find(beh==1 &  Position==4);
    min_trials=min([length(a1) length(a2) length(a3) length(a4)]);
    a1 = a1(randperm(length(a1)));
    a2 = a2(randperm(length(a2)));
    a3 = a3(randperm(length(a3)));
    a4 = a4(randperm(length(a4)));
    a1=a1(1:min_trials);
    a2=a2(1:min_trials);
    a3=a3(1:min_trials);
    a4=a4(1:min_trials);
    POS{1,1}=a1;
    POS{1,2}=a2;
    POS{1,3}=a3;
    POS{1,4}=a4;
    winpre=100:350; %% Pre-stimulus period
    winpost=400:800; %% Post stimulus period
    PVAL=[];
    pval_all=[];
    X_all=[];
    for channel=1:48;
        k_chan=k_chan+1;
        SIG_cue=[];
        SIG_target=[];
        for pos=1:4;
            tmp=POS{1,pos};
            A=squeeze(mean(dataC(channel,winpre,tmp),2));
            B1=squeeze(mean(dataC(channel,winpost,tmp),2));
            [p,h]=ranksum(A,B1);
            SIG_cue=[SIG_cue p];
            A=squeeze(mean(dataT(channel,winpre,tmp),2));
            B2=squeeze(mean(dataT(channel,winpost,tmp),2));
            [p,h]=ranksum(A,B2);
            SIG_target=[SIG_target p];
        end
        if min(SIG_cue)<0.05 & min(SIG_target)<0.05;
            sel_chan=sel_chan+1;
            selective_channel(sel_chan) = k_chan;
            [SIG_cue cue]=sort(SIG_cue,'ascend');
            [SIG_target target]=sort(SIG_target,'ascend');
            A=squeeze(mean(dataC(channel,winpost,POS{1,target(1)}),2));%% Prefered
            B=squeeze(mean(dataC(channel,winpost,POS{1,target(4)}),2));%% No Prefered
            p = ranksum(A,B);
            IND(k_chan)= (mean(A)-mean(B))./(mean(A)+mean(B));
            SIG_tot(k_chan)=p;
        end
    end
    clear data*
end

%%% PLOT Modulation index

I=IND(selective_channel);
P=SIG_tot(selective_channel);
w=min(I):(max(I)-min(I))/125:max(I);

clear COUNT*
k=0;
for j=2:length(w)-1;
    k=k+1;
    COUNT(k)=length(find(I<w(j+1) & I>w(j-1)))/length(I);
    a=find(I<w(j+1) & I>w(j-1));
    PP=P(a);
    COUNT_sig(k)=length(find(PP<0.05 & PP>0))/length(I);
end
%figure, bar(w(2:end-1),COUNT,'k'), hold all,  bar(w(2:end-1),COUNT_sig,'w')
figure; plot(w(2:end-1),COUNT_sig,'b'); hold all; plot(w(2:end-1),COUNT,'r');xlabel('Modulation Index'); ylabel('Proportion of cells'); legend('Significant (p<0.05)','No significant')


%%%%%% 

%%%FIGURE 2C MUA PREFERED - LEAST-PREFERED


file = dir('*spikes.mat');
kkk = 0
for session=[1:length(file)]
    kkk=kkk+1;
    load(file(session).name);
    
    sm = 60; % Smoothing parameter
    window = ones(sm,1)/sm;
    time=-1500:500;
    winT=find(time==-500):find(time==500);
    time=-500:1500;
    winC=find(time==-500):find(time==500);
    for channel=1:48;
        for trial=1:size(dataT,3);
            tr=conv(squeeze(dataC(channel,winC,trial)),window,'same');
            tr=tr(sm/2:end-sm/2);
            % a=resample(a,100,1000);
            % dataC2(channel,:,trial)=(a-mean(a(1:300)))./std(a(1:300));
            dataC_smooth(channel,:,trial)=tr;
            tr=conv(squeeze(dataT(channel,winT,trial)),window,'same');
            tr=tr(sm/2:end-sm/2);
            % a=resample(a,100,1000);
            %dataT2(channel,:,trial)=(a1-mean(a(1:300)))./std(a(1:300));;
            dataT_smooth(channel,:,trial)=tr;
        end
    end
    beh = Behavior(:,1);  % 1:hit, 2:miss,
    Position = Behavior(:,2);
    a1 = find(beh==1 &  Position==1);
    a2 = find(beh==1 &  Position==2);
    a3 = find(beh==1 &  Position==3);
    a4 = find(beh==1 &  Position==4);
    min_trials=min([length(a1) length(a2) length(a3) length(a4)]);
    a1 = a1(randperm(length(a1)));
    a2 = a2(randperm(length(a2)));
    a3 = a3(randperm(length(a3)));
    a4 = a4(randperm(length(a4)));
    a1=a1(1:min_trials);
    a2=a2(1:min_trials);
    a3=a3(1:min_trials);
    a4=a4(1:min_trials);
    POS{1,1}=a1;
    POS{1,2}=a2;
    POS{1,3}=a3;
    POS{1,4}=a4;
    winpre=100:350; %% Pre-stimulus period
    winpost=400:800; %% Post stimulus period
    PVAL=[];
    pval_all=[];
    X_all=[];
    for channel=1:48
        SFR=[];
        for pos=1:4
            tmp=POS{1,pos};
            A=squeeze(mean(dataC_smooth(channel,winpre,tmp),2));
            B=squeeze(mean(dataC_smooth(channel,winpost,tmp),2));
            [p,h]=ranksum(A,B);
            SFR=[SFR p];
        end
        [SFRs,SFRind]=sort(SFR);
        PREF_cue(channel,kkk,:)=(squeeze(mean(dataC_smooth(channel,:,POS{1,SFRind(1)}),3))-squeeze(mean(mean(dataC_smooth(channel,winpre,POS{1,SFRind(1)}),3),2)))./std(squeeze(mean(dataC_smooth(channel,winpre,POS{1,SFRind(1)}),3)));
        NO_PREF_cue(channel,kkk,:)=(squeeze(mean(dataC_smooth(channel,:,POS{1,SFRind(4)}),3))-squeeze(mean(mean(dataC_smooth(channel,winpre,POS{1,SFRind(4)}),3),2)))./std(squeeze(mean(dataC_smooth(channel,winpre,POS{1,SFRind(4)}),3)));
        PREF_target(channel,kkk,:)=(squeeze(mean(dataT_smooth(channel,:,POS{1,SFRind(1)}),3))-squeeze(mean(mean(dataC_smooth(channel,winpre,POS{1,SFRind(1)}),3),2)))./std(squeeze(mean(dataC_smooth(channel,winpre,POS{1,SFRind(1)}),3)));
        NO_PREF_target(channel,kkk,:)=(squeeze(mean(dataT_smooth(channel,:,POS{1,SFRind(4)}),3))-squeeze(mean(mean(dataC_smooth(channel,winpre,POS{1,SFRind(4)}),3),2)))./std(squeeze(mean(dataC_smooth(channel,winpre,POS{1,SFRind(4)}),3)));
        
    end
    clear data*
end

PREF_target2 = reshape(PREF_target,48*18,942);
PREF_cue2 = reshape(PREF_cue,48*18,942);
NO_PREF_target2 = reshape(NO_PREF_target,48*18,942);
NO_PREF_cue2 = reshape(NO_PREF_cue,48*18,942);

figure, semshade(PREF_cue2,0.1,[0 0 0],1:942,10,0,1); hold all, semshade(NO_PREF_cue2,0.1,[0.6 0.6 0.6],1:942,10,0,1); 

figure, semshade(PREF_target2,0.1,[0 0 0],1:942,10,0,1); hold all, semshade(NO_PREF_target2,0.1,[0.6 0.6 0.6],1:942,10,0,1);



    
   


