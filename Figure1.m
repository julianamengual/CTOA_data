path_data = cd('../Behavior/');
file = dir([path_data '/', '*.mat']);

for session=1:length(file);
    
    load(file(session).name,'ARRAY')
    %% Array is a matrix N x C, N: Number of trials, C = Behavioral Parameters
    
    % 1: Trial Number
    % 2: x-coordinate of the cue (1, -1)
    % 3: y-coordinate of the cue (1, -1)
    % 4: 1: There is distractor in the trial, 0: There is no distractor in the trial
    % 5: Correct trial (the monkey responded to the target)
    % 7: Time passed between Cue and distractor
    % 8: Time passed between Distractor and Target
    % 9: Time passed between the response and the distractor
    %10: Time passed between the response and the Cue
    %11: Time passed between the response and the target (The reaction
    %time)
    %12: Time between the Cue and the Target (CTOA)
    %13: Behavior: 1:hit, 2:miss, 3: False alarm, 4: Abort trial
    
    hits=find(ARRAY(:,13)==1 );
    hits_and_misses=find((ARRAY(:,13)==1 | ARRAY(:,13)==2));
    hit_rate(session)=length(hits)/length(hits_and_misses); % Hit rate
    RT_hits(session) = mean(ARRAY(hits,11)); %% Reaction time in hits
    

    
    false_alarms=find(ARRAY(:,13)==3 & ARRAY(:,9)>200 & ARRAY(:,9)<900 & ARRAY(:,4)>0);
    trials_with_distractors=find(ARRAY(:,4)>0);
    RT_false_alarms(session) = mean(ARRAY(false_alarms,9)); % Reaction time in false alarms
    false_alarm_rate(session)=length(false_alarms)/length(trials_with_distractors); % False alarm rate
   
    
end





