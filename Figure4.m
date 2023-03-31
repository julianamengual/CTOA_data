
d=dir([cd '/' '*data.mat']);
kkk=0;
k_chan=0;
sel_chan=0;


for iii=[1:length(d)]
    iii
    load(d(iii).name),
    
    tmp=find(Behavior(:,1)==1);
    Behavior=Behavior(tmp,:);
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    
    
    [mu1, sig1, betas1,beta_fast1,beta_main1, graphhandles1, x_1,y_1] = reciprobit2(Behavior(tmp,4)/1000);
    [res_x, idx_of_result1] = knee_pt(x_1,y_1,0);
    a11=find(Behavior(:,1)==1 & Behavior(:,4)<-1000*(1./x_1(idx_of_result1)) );
    Behavior(a11,:)=[];
    dataT(:,:,a11)=[];
    dataC(:,:,a11)=[];
    
    Behavior(:,9)=Behavior(:,9).*14*sqrt(2)/2;
    a1 = find(Behavior(:,1)==1 &  (Behavior(:,2)==1));
    a2 = find(Behavior(:,1)==1 &  (Behavior(:,2)==2));
    a3 = find(Behavior(:,1)==1 &  (Behavior(:,2)==3));
    a4 = find(Behavior(:,1)==1 &  (Behavior(:,2)==4));
    a=min([length(a1) length(a2) length(a3) length(a4)]);
    
    b=randperm(a);
    a1=a1(b);a2=a2(b);a3=a3(b);a4=a4(b);
    a1=a1(1:a);
    a2=a2(1:a);
    a3=a3(1:a);
    a4=a4(1:a);
    TMP{1,1}=a1;
    TMP{1,2}=a2;
    TMP{1,3}=a3;
    TMP{1,4}=a4;
    winpre=300:400;
    winpost=400:800;
    PVAL=[];
    pval_all=[];
    X_all=[];
    for channel=1:48;
        k_chan=k_chan+1;
        SIG_cue=[];
        SIG_target=[];
        for pos=1:4;
            tmp=TMP{1,pos};
            A=squeeze(mean(dataC(channel,winpre,tmp),2));
            B1=squeeze(mean(dataC(channel,winpost,tmp),2));
            [h, p]=ttest(A,B1);
            SIG_cue=[SIG_cue p];
            A=squeeze(mean(dataT(channel,winpre,tmp),2));
            B2=squeeze(mean(dataT(channel,winpost,tmp),2));
            [h, p]=ttest(A,B2);
            SIG_target=[SIG_target p];
        end
        if min(SIG_cue)<0.05  & min(SIG_target)<0.05;
            sel_chan=sel_chan+1
            
            
            selective_channel2(sel_chan)=k_chan;
            [~,pref_pos] = min(SIG_target);
            clear A
            array_trials = 1:size(Behavior,1);
            
            clear A
            k=0;
            for t=1:10:801-10
                k=k+1;
                for ii=1:10
                    trial{ii}=find(Behavior(:,3)>prctile(Behavior(:,3),10*(ii-1)) & Behavior(:,3)<prctile(Behavior(:,3),10*ii));
                    
                end
                for i=1:length(trial);
                    num_trial(i)=length(trial{i});
                end;
                min_num_trial=min(num_trial);
                
                for ii=1:10
                    x=trial{ii};
                    A(:,ii) = squeeze(mean(dataT(channel,t:t+10,x(1:min_num_trial)),2));
                end
                [p]=friedman(A);
                
                clear A
                SIG_CTOA_pref(sel_chan,k)=p;
                
                
            end
            
            k=0;
            for t=1:10:801-10
                k=k+1;
                for ii=1:10
                    trial{ii}=find(Behavior(:,4)>prctile(Behavior(:,4),10*(ii-1)) & Behavior(:,4)<prctile(Behavior(:,4),10*ii));
                    
                end
                for i=1:length(trial);
                    num_trial(i)=length(trial{i});
                end;
                min_num_trial=min(num_trial);
                for ii=1:10
                    x=trial{ii};
                    A(:,ii) = squeeze(mean(dataT(channel,t:t+10,x(1:min_num_trial)),2));
                end
                [p]=friedman(A);
                
                clear A
                SIG_RT_pref(sel_chan,k)=p;
            end
            k=0;
            for t=1:10:801-10
                k=k+1;
                for ii=1:10
                    trial{ii}=find(Behavior(:,9)>prctile(Behavior(:,9),10*(ii-1)) & Behavior(:,9)<prctile(Behavior(:,9),10*ii));
                    
                end
                for i=1:length(trial);
                    num_trial(i)=length(trial{i});
                end;
                min_num_trial=min(num_trial);
                for ii=1:10
                    x=trial{ii};
                    A(:,ii) = squeeze(mean(dataT(channel,t:t+10,x(1:min_num_trial)),2));
                end
                [p]=friedman(A);
                
                clear A
                SIG_TA_pref(sel_chan,k)=p;
            end
            
            k=0;
            for t=1:10:801-10
                k=k+1;
                a1 = find(Behavior(:,1)==1 & Behavior(:,2)==1);
                a2 = find(Behavior(:,1)==1 & Behavior(:,2)==2);
                a3 = find(Behavior(:,1)==1 & Behavior(:,2)==3);
                a4 = find(Behavior(:,1)==1 & Behavior(:,2)==4);
                a=min([length(a1) length(a2) length(a3) length(a4)]);
                b=randperm(a);
                a1=a1(b);a2=a2(b);a3=a3(b);a4=a4(b);
                
                POS(:,1)=squeeze(mean(dataT(channel,t:t+10,a1),2));
                POS(:,2)=squeeze(mean(dataT(channel,t:t+10,a2),2));
                POS(:,3)=squeeze(mean(dataT(channel,t:t+10,a3),2));
                POS(:,4)=squeeze(mean(dataT(channel,t:t+10,a4),2));
                SIG_POS(sel_chan,k)=friedman(POS);
                
            end
            clear POS
        end
    end
    clear data*
    close all
    
end

for t=1:80, 
    sel_ctoa(t) = length(find(SIG_CTOA_pref(:,t)<0.05))/size(SIG_CTOA_pref,1);
    sel_ta(t) = length(find(SIG_TA_pref(:,t)<0.05))/size(SIG_TA_pref,1);
    sel_rt(t) = length(find(SIG_RT_pref(:,t)<0.05))/size(SIG_RT_pref,1);
    sel_position(t) = length(find(SIG_POS(:,t)<0.05))/size(SIG_POS,1);
end
    
figure, plot(sel_ctoa);hold all,  plot(sel_ta); hold all,  plot(sel_rt); hold all, plot(sel_position); 
legend('CTOA','TA','RT','Position')
    
    