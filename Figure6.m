clear all
d=dir([cd '/' '*data.mat']);
kkk=0;
%firingRates=zeros(18,48,3,3,3,801,57);
for iii=[1:length(d)]
    kkk=kkk+1
    load(d(iii).name,'data*','Behavior');
    iii
    tmp=find(Behavior(:,1)==1 & Behavior(:,3)<3500 ) ;
    % We categorize the CTOA in 4 different ranges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ==== Smoothing ==== %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Behavior=Behavior(tmp,:);
    dataT=dataT(:,:,tmp);
    dataC=dataC(:,:,tmp);
    Behavior(:,9)=Behavior(:,9).*14*sqrt(2)/2;
    a1=find(Behavior(:,3)<prctile(Behavior(:,3),33));
    a2=find(Behavior(:,3)<prctile(Behavior(:,3),66) & Behavior(:,3)>prctile(Behavior(:,3),33) );
    a3=find(Behavior(:,3)>prctile(Behavior(:,3),66));
    Behavior_early=Behavior(a1,:);
    Behavior_middle=Behavior(a2,:);
    Behavior_late=Behavior(a3,:);
    dataT_early=dataT(:,1:2:end,a1);
    dataT_middle=dataT(:,1:2:end,a2);
    dataT_late=dataT(:,1:2:end,a3);
    
    a11=find(Behavior_early(:,9)<prctile(Behavior_early(:,9),33));
    a13=find(Behavior_early(:,9)>prctile(Behavior_early(:,9),66));
    a12=find(Behavior_early(:,9)<prctile(Behavior_early(:,9),66) & Behavior_early(:,9)>prctile(Behavior_early(:,9),33) );
    
    a31=find(Behavior_late(:,9)<prctile(Behavior_late(:,9),33));
    a33=find(Behavior_late(:,9)>prctile(Behavior_late(:,9),66));
    a32=find(Behavior_late(:,9)<prctile(Behavior_late(:,9),66) & Behavior_late(:,9)>prctile(Behavior_late(:,9),33) );
    
    a21=find(Behavior_middle(:,9)<prctile(Behavior_middle(:,9),33));
    a23=find(Behavior_middle(:,9)>prctile(Behavior_middle(:,9),66));
    a22=find(Behavior_middle(:,9)<prctile(Behavior_middle(:,9),66) & Behavior_middle(:,9)>prctile(Behavior_middle(:,9),33) );
    
    Behavior_early_short=Behavior_early(a11,:);
    Behavior_early_medium=Behavior_early(a12,:);
    Behavior_early_long=Behavior_early(a13,:);
    
    
    Behavior_middle_short=Behavior_middle(a21,:);
    Behavior_middle_medium=Behavior_middle(a22,:);
    Behavior_middle_long=Behavior_middle(a23,:);
    
    Behavior_late_short=Behavior_late(a31,:);
    Behavior_late_medium=Behavior_late(a32,:);
    Behavior_late_long=Behavior_late(a33,:);
    
    
    dataT_early_short=dataT_early(:,:,a11);
    dataT_early_medium=dataT_early(:,:,a12);
    dataT_early_long=dataT_early(:,:,a13);
    
    dataT_middle_short=dataT_middle(:,:,a21);
    dataT_middle_medium=dataT_middle(:,:,a22);
    dataT_middle_long=dataT_middle(:,:,a23);
    
    dataT_late_short=dataT_late(:,:,a31);
    dataT_late_medium=dataT_late(:,:,a32);
    dataT_late_long=dataT_late(:,:,a33);
    
    
    a111=find(Behavior_early_short(:,4)<prctile(Behavior_early_short(:,4),33));
    a112=find(Behavior_early_short(:,4)>prctile(Behavior_early_short(:,4),33) &Behavior_early_short(:,4)<prctile(Behavior_early_short(:,4),66));
    a113=find(Behavior_early_short(:,4)>prctile(Behavior_early_short(:,4),66));
    
    a121=find(Behavior_early_medium(:,4)<prctile(Behavior_early_medium(:,4),33));
    a122=find(Behavior_early_medium(:,4)>prctile(Behavior_early_medium(:,4),33) &Behavior_early_medium(:,4)<prctile(Behavior_early_medium(:,4),66));
    a123=find(Behavior_early_medium(:,4)>prctile(Behavior_early_medium(:,4),66));
    
    a131=find(Behavior_early_long(:,4)<prctile(Behavior_early_long(:,4),33));
    a132=find(Behavior_early_long(:,4)>prctile(Behavior_early_long(:,4),33) &Behavior_early_long(:,4)<prctile(Behavior_early_long(:,4),66));
    a133=find(Behavior_early_long(:,4)>prctile(Behavior_early_long(:,4),66));
    
    a211=find(Behavior_middle_short(:,4)<prctile(Behavior_middle_short(:,4),33));
    a212=find(Behavior_middle_short(:,4)>prctile(Behavior_middle_short(:,4),33) &Behavior_middle_short(:,4)<prctile(Behavior_middle_short(:,4),66));
    a213=find(Behavior_middle_short(:,4)>prctile(Behavior_middle_short(:,4),66));
    
    a221=find(Behavior_middle_medium(:,4)<prctile(Behavior_middle_medium(:,4),33));
    a222=find(Behavior_middle_medium(:,4)>prctile(Behavior_middle_medium(:,4),33) &Behavior_middle_medium(:,4)<prctile(Behavior_middle_medium(:,4),66));
    a223=find(Behavior_middle_medium(:,4)>prctile(Behavior_middle_medium(:,4),66));
    
    a231=find(Behavior_middle_long(:,4)<prctile(Behavior_middle_long(:,4),33));
    a232=find(Behavior_middle_long(:,4)>prctile(Behavior_middle_long(:,4),33) &Behavior_middle_long(:,4)<prctile(Behavior_middle_long(:,4),66));
    a233=find(Behavior_middle_long(:,4)>prctile(Behavior_middle_long(:,4),66));
    
    a311=find(Behavior_late_short(:,4)<prctile(Behavior_late_short(:,4),33));
    a312=find(Behavior_late_short(:,4)>prctile(Behavior_late_short(:,4),33) &Behavior_late_short(:,4)<prctile(Behavior_late_short(:,4),66));
    a313=find(Behavior_late_short(:,4)>prctile(Behavior_late_short(:,4),66));
    
    a321=find(Behavior_late_medium(:,4)<prctile(Behavior_late_medium(:,4),33));
    a322=find(Behavior_late_medium(:,4)>prctile(Behavior_late_medium(:,4),33) &Behavior_late_medium(:,4)<prctile(Behavior_late_medium(:,4),66));
    a323=find(Behavior_late_medium(:,4)>prctile(Behavior_late_medium(:,4),66));
    
    a331=find(Behavior_late_long(:,4)<prctile(Behavior_late_long(:,4),33));
    a332=find(Behavior_late_long(:,4)>prctile(Behavior_late_long(:,4),33) &Behavior_late_long(:,4)<prctile(Behavior_late_long(:,4),66));
    a333=find(Behavior_late_long(:,4)>prctile(Behavior_late_long(:,4),66));
    
    dataT_early_short_fast=dataT_early_short(:,:,a111);
    dataT_early_short_nofast=dataT_early_short(:,:,a112);
    dataT_early_short_slow=dataT_early_short(:,:,a113);
    
    Matr_condT(kkk,:,1,1,1,:)=mean(dataT_early_short_fast,3);
    firingRates(kkk,:,1,1,1,:,1:length(a111))=dataT_early_short_fast;
    Matr_condT(kkk,:,1,1,2,:)=mean(dataT_early_short_nofast,3);
    firingRates(kkk,:,1,1,2,:,1:length(a112))=dataT_early_short_nofast;
    Matr_condT(kkk,:,1,1,3,:)=mean(dataT_early_short_slow,3);
    firingRates(kkk,:,1,1,3,:,1:length(a113))=dataT_early_short_slow;
    
    dataT_early_medium_fast=dataT_early_medium(:,:,a121);
    dataT_early_medium_nofast=dataT_early_medium(:,:,a122);
    dataT_early_medium_slow=dataT_early_medium(:,:,a123);
    
    
    Matr_condT(kkk,:,1,2,1,:)=mean(dataT_early_medium_fast,3);
    firingRates(kkk,:,1,2,1,:,1:length(a121))=dataT_early_medium(:,:,a121);
    Matr_condT(kkk,:,1,2,2,:)=mean(dataT_early_medium_nofast,3);
    firingRates(kkk,:,1,2,2,:,1:length(a122))=dataT_early_medium(:,:,a122);
    Matr_condT(kkk,:,1,2,3,:)=mean(dataT_early_medium_slow,3);
    firingRates(kkk,:,1,2,3,:,1:length(a123))=dataT_early_medium(:,:,a123);
    
    
    dataT_early_long_fast=dataT_early_long(:,:,a131);
    dataT_early_long_nofast=dataT_early_long(:,:,a132);
    dataT_early_long_slow=dataT_early_long(:,:,a133);
    
    
    Matr_condT(kkk,:,1,3,1,:)=mean(dataT_early_long_fast,3);
    firingRates(kkk,:,1,3,1,:,1:length(a131))=dataT_early_long(:,:,a131);
    Matr_condT(kkk,:,1,3,2,:)=mean(dataT_early_long_nofast,3);
    firingRates(kkk,:,1,3,2,:,1:length(a132))=dataT_early_long(:,:,a132);
    Matr_condT(kkk,:,1,3,3,:)=mean(dataT_early_long_slow,3);
    firingRates(kkk,:,1,3,3,:,1:length(a133))=dataT_early_long(:,:,a133);
    
    dataT_middle_short_fast=dataT_middle_short(:,:,a211);
    dataT_middle_short_nofast=dataT_middle_short(:,:,a212);
    dataT_middle_short_slow=dataT_middle_short(:,:,a213);
    
    Matr_condT(kkk,:,2,1,1,:)=mean(dataT_middle_short_fast,3);
    firingRates(kkk,:,2,1,1,:,1:length(a211))=dataT_middle_short(:,:,a211);
    Matr_condT(kkk,:,2,1,2,:)=mean(dataT_middle_short_nofast,3);
    firingRates(kkk,:,2,1,2,:,1:length(a212))=dataT_middle_short(:,:,a212);
    Matr_condT(kkk,:,2,1,3,:)=mean(dataT_middle_short_slow,3);
    firingRates(kkk,:,2,1,3,:,1:length(a213))=dataT_middle_short(:,:,a213);
    
    dataT_middle_medium_fast=dataT_middle_medium(:,:,a221);
    dataT_middle_medium_nofast=dataT_middle_medium(:,:,a222);
    dataT_middle_medium_slow=dataT_middle_medium(:,:,a223);
    
    Matr_condT(kkk,:,2,2,1,:)=mean(dataT_middle_medium_fast,3);
    firingRates(kkk,:,2,2,1,:,1:length(a221))=dataT_middle_medium(:,:,a221);
    Matr_condT(kkk,:,2,2,2,:)=mean(dataT_middle_medium_nofast,3);
    firingRates(kkk,:,2,2,2,:,1:length(a222))=dataT_middle_medium(:,:,a222);
    Matr_condT(kkk,:,2,2,3,:)=mean(dataT_middle_medium_slow,3);
    firingRates(kkk,:,2,2,3,:,1:length(a223))=dataT_middle_medium(:,:,a223);
    
    dataT_middle_long_fast=dataT_middle_long(:,:,a231);
    dataT_middle_long_nofast=dataT_middle_long(:,:,a232);
    dataT_middle_long_slow=dataT_middle_long(:,:,a233);
    
    Matr_condT(kkk,:,2,3,1,:)=mean(dataT_middle_long_fast,3);
    firingRates(kkk,:,2,3,1,:,1:length(a231))=dataT_middle_long(:,:,a231);
    Matr_condT(kkk,:,2,3,2,:)=mean(dataT_middle_long_nofast,3);
    firingRates(kkk,:,2,3,2,:,1:length(a232))=dataT_middle_long(:,:,a232);
    Matr_condT(kkk,:,2,3,3,:)=mean(dataT_middle_long_slow,3);
    firingRates(kkk,:,2,3,3,:,1:length(a233))=dataT_middle_long(:,:,a233);
    
    dataT_late_short_fast=dataT_late_short(:,:,a311);
    dataT_late_short_nofast=dataT_late_short(:,:,a312);
    dataT_late_short_slow=dataT_late_short(:,:,a313);
    
    Matr_condT(kkk,:,3,1,1,:)=mean(dataT_late_short_fast,3);
    firingRates(kkk,:,3,1,1,:,1:length(a311))=dataT_late_short(:,:,a311);
    Matr_condT(kkk,:,3,1,2,:)=mean(dataT_late_short_nofast,3);
    firingRates(kkk,:,3,1,2,:,1:length(a312))=dataT_late_short(:,:,a312);
    Matr_condT(kkk,:,3,1,3,:)=mean(dataT_late_short_slow,3);
    firingRates(kkk,:,3,1,3,:,1:length(a313))=dataT_late_short(:,:,a313);
    
    dataT_late_medium_fast=dataT_late_medium(:,:,a321);
    dataT_late_medium_nofast=dataT_late_medium(:,:,a322);
    dataT_late_medium_slow=dataT_late_medium(:,:,a323);
    
    Matr_condT(kkk,:,3,2,1,:)=mean(dataT_late_medium_fast,3);
    firingRates(kkk,:,3,2,1,:,1:length(a321))=dataT_late_medium(:,:,a321);
    Matr_condT(kkk,:,3,2,2,:)=mean(dataT_late_medium_nofast,3);
    firingRates(kkk,:,3,2,2,:,1:length(a322))=dataT_late_medium(:,:,a322);
    Matr_condT(kkk,:,3,2,3,:)=mean(dataT_late_medium_slow,3);
    firingRates(kkk,:,3,2,3,:,1:length(a323))=dataT_late_medium(:,:,a323);
    
    dataT_late_long_fast=dataT_late_long(:,:,a331);
    dataT_late_long_nofast=dataT_late_long(:,:,a332);
    dataT_late_long_slow=dataT_late_long(:,:,a333);
    
    Matr_condT(kkk,:,3,3,1,:)=mean(dataT_late_long_fast,3);
    firingRates(kkk,:,3,3,1,:,1:length(a331))=dataT_late_long(:,:,a331);
    Matr_condT(kkk,:,3,3,2,:)=mean(dataT_late_long_nofast,3);
    firingRates(kkk,:,3,3,2,:,1:length(a332))=dataT_late_long(:,:,a332);
    Matr_condT(kkk,:,3,3,3,:)=mean(dataT_late_long_slow,3);
    firingRates(kkk,:,3,3,3,:,1:length(a333))=dataT_late_long(:,:,a333);
    
    
    
    clear dataT* Behavior*
    
    
    
    clear data*
end


Matr_CondT=reshape(Matr_condT,18*48,3,3,3,401);
firingRates=reshape(firingRates,18*48,3,3,3,401,size(firingRates,7));
combinedParams = {{1, [1 4]}, {2, [2 4]},{3, [3 4]}, {4}};

margNames = {'CTOA', 'TA', 'RT' ,'Condition-independent'};
margColours = jet(4);

DATA2=Matr_CondT(:,[1 3],[1:3],[1 3],1:400);
DATA2_all=firingRates(:,[1 3],[1:3],[1 3],1:400,:);

trialNum=size(firingRates,6);
X = DATA2(:,:);
X = bsxfun(@minus, X, mean(X,2));
[W,s,~] = svd(X, 'econ');
W = W(:,1:48);

[W,V,whichMarg] = dpca(DATA2, 864, ...
    'combinedParams', combinedParams);
explVar = dpca_explainedVariance(DATA2, W, V, ...
    'combinedParams', combinedParams);
dpca_plot(DATA2, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', [],                        ...
    'timeEvents', [],               ...
    'timeMarginalization',3, ...
    'legendSubplot', 16);
optimalLambda = dpca_optimizeLambda(DATA2(:,:,:,:,1:40:end), DATA2_all(:,:,:,:,1:40:end,:), size(firingRates,6), ...
    'combinedParams', combinedParams, ...
    'simultaneous', true, ...
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'filename', []);

Cnoise = dpca_getNoiseCovariance(DATA2, ...
    DATA2_all, trialNum, 'simultaneous', true);
[W,V,whichMarg] = dpca(DATA2, 864, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(DATA2, W, V, ...
    'combinedParams', combinedParams);
dpca_plot(DATA2, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', [],                        ...
    'timeEvents', [],               ...
    'timeMarginalization',3, ...
    'legendSubplot', 16);

S=2
decodingClasses = {[(1:S)']};
combinedParams = {{1, [1 2]}, {2}};

W_target = W;

C=[1 3 4 2]
figure,
for kkk=1:length(C)
    subplot(1,4,kkk);
    comp2=C(kkk)
    P=X'*W_target(:,comp2);
    PP=reshape(P,1,2,3,2,400);
    
    
    
    comp=1
    
    
    for i=[1 2]
        for j=1:3
            for k=1:2
                p=plot(squeeze(PP(comp,i,j,k,1:400)));
                p.LineWidth=2
                if k==1
                    p.LineStyle='--';
                    if i==1 & j==1;
                        p.Color = [136 8 0]./255;
                    end
                    if i==1 & j==2
                        p.Color = [255 0 0]./255;
                    end
                    if i==1 & j==3
                        
                        p.Color = [255, 203, 192]./255;
                    end
                    
                    if i==2 & j==1;
                        p.Color = [46,139,87]./255;
                    end
                    if i==2 & j==2
                        p.Color = [50,205,50]./255;
                    end
                    
                    if i==2 & j==3
                        
                        p.Color = [124,252,0]./255;
                    end
                    hold on;
                else
                    if i==1 & j==1;
                        p.Color = [136 8 0]./255;
                    end
                    if i==1 & j==2
                        p.Color = [255 0 0]./255;
                    end
                    if i==1 & j==3
                        
                        p.Color = [255, 203, 192]./255;
                    end
                    
                    if i==2 & j==1;
                        p.Color = [46,139,87]./255;
                    end
                    if i==2 & j==2
                        p.Color = [50,205,50]./255;
                    end
                    
                    if i==2 & j==3
                        
                        p.Color = [124,252,0]./255;
                    end
                    hold on,
                end
                
                
            end
        end
    end
end




subplot(1,2,2);
for i=[1 3]
    for j=1:3
        for k=1:2
            p=plot(squeeze(PP(comp,i,j,k,82:end)));
            p.LineWidth=2
            if k==1
                p.LineStyle='--';
                if i==1 & j==1;
                    p.Color = [136 8 0]./255;
                end
                if i==1 & j==2
                    p.Color = [255 0 0]./255;
                end
                if i==1 & j==3
                    
                    p.Color = [255, 203, 192]./255;
                end
                
                if i==3 & j==1;
                    p.Color = [46,139,87]./255;
                end
                if i==3 & j==2
                    p.Color = [50,205,50]./255;
                end
                
                if i==3 & j==3
                    
                    p.Color = [124,252,0]./255;
                end
                hold on;
            else
                if i==1 & j==1;
                    p.Color = [136 8 0]./255;
                end
                if i==1 & j==2
                    p.Color = [255 0 0]./255;
                end
                if i==1 & j==3
                    
                    p.Color = [255, 203, 192]./255;
                end
                
                if i==3 & j==1;
                    p.Color = [46,139,87]./255;
                end
                if i==3 & j==2
                    p.Color = [50,205,50]./255;
                end
                
                if i==3 & j==3
                    
                    p.Color = [124,252,0]./255;
                end
                hold on,
            end
            
            
        end
    end
end
axis([ 1 80 -0.2 0.2])








