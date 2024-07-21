function [] = filterAP2neg_3ch(data, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', @isstruct);
ip.addParameter('mask',[],@isnumeric);
ip.addParameter('AP2ch',2,@isnumeric);
ip.addParameter('CLTAch',3,@isnumeric);
ip.addParameter('intensityThreshold',0 ,@isnumeric);
ip.addParameter('SNR',[1.05 1.05],@isnumeric);
ip.addParameter('A_SNR',[1.08 1.3],@isnumeric); %Artificial definition SNR
ip.addParameter('Lifetime_threshold',[5 Inf],@isnumeric);%minlifetime 5s maxlifetime 120s
ip.addParameter('GAP_ratio',0.2,@isnumeric); % permit 0.2 GAP
ip.parse(varargin{:});
snr = ip.Results.SNR;
asnr = ip.Results.A_SNR;
gapratio = ip.Results.GAP_ratio;
intensityThreshold = ip.Results.intensityThreshold;
lifetime_t = ip.Results.Lifetime_threshold;
PLC_AP2neg_summary = struct;
PLC_AP2negLong_summary = struct;
PLC_AP2pos_summary = struct;
cellArea = [];
parfor nd = 1:length(data)
    fprintf('filter start %d movies\n',nd);
    temtracks = load([data(nd).source 'Tracking\ProcessedTracks.mat']);
    tracks = temtracks.tracks;
    [mask,cellArea(nd)] = load_cellmask(data(nd),ip.Results.mask);   
    % sort lifetime for cmedataviewer2_KM
    [~,Index] = sort([tracks.lifetime_s], 'descend');
    AP2neg_idx=[];
    AP2neg_idx_long =[];
    AP2pos_idx=[];
    labels = zeros(1,length(tracks));
    % prefilter
    idx = [];
    for i = 1:length(tracks)
        Signalmch = tracks(i).A(1,:)+tracks(i).c(1,:);
        Signalmchlow = tracks(i).A(1,:)+tracks(i).c(1,:)-tracks(i).A_pstd(1,:);
        Noisehigh = tracks(i).c(1,:)+1.96*tracks(i).sigma_r(1,:);
        Noisemch = tracks(i).c(1,:);
        real_detectionframe = find(tracks(i).gapVect==0);
        maxSignalmch = max(tracks(i).A(1,real_detectionframe)+tracks(i).c(1,real_detectionframe));
        tracks_meanSNR = mean(Signalmch)/mean(Noisemch);
        tracks_maxSNR = maxSignalmch/mean(Noisemch);
        x_position = round(nanmean(tracks(i).x(1,:)));
        y_position = round(nanmean(tracks(i).y(1,:)));
        gapnumber = length(find(tracks(i).gapVect==1));
        maskfilter = mask(y_position, x_position) == 1;
        snrfilter = tracks_meanSNR>=snr(1) && tracks_maxSNR>=snr(2);
        a_defined_filter = mean(Signalmchlow)/mean(Noisehigh) >asnr(1);
        max_firstfilter = maxSignalmch/(tracks(i).A(1,1)+tracks(i).c(1,1)) > asnr(2);
        gapfilter = gapnumber/length(tracks(i).t) < gapratio;
        lifetimefilter = tracks(i).lifetime_s>=lifetime_t(1) && tracks(i).lifetime_s<=lifetime_t(2);
        intensityfilter = max(tracks(i).A(1,:))>=intensityThreshold;
        if  maskfilter && snrfilter && gapfilter && a_defined_filter && max_firstfilter && lifetimefilter && intensityfilter
            idx = [idx i];
        end
    end
    tracks = tracks(idx);
    [tmp_neg_index,tmp_neg_longindex,tmp_pos_index]=deal([]); 
    for i = 1:length(tracks)
        if tracks(i).significantSlave(ip.Results.AP2ch)==0 ...,
                && tracks(i).significantSlave(ip.Results.CLTAch)==1 ...,
                && (tracks(i).catIdx ==1 || tracks(i).catIdx ==2)
            AP2neg_idx = [AP2neg_idx idx(i)];
            tmp_neg_index = [tmp_neg_index i];
        end
        if tracks(i).significantSlave(ip.Results.AP2ch)==0 ...,
                && tracks(i).significantSlave(ip.Results.CLTAch)==1 ...,
                && (tracks(i).catIdx ==3 || tracks(i).catIdx ==4)
            AP2neg_idx_long=[AP2neg_idx_long idx(i)];
            tmp_neg_longindex = [tmp_neg_longindex i];
        end
        if tracks(i).significantSlave(ip.Results.AP2ch)==1 ...,
                && tracks(i).significantSlave(ip.Results.CLTAch)==1 ...,
                && (tracks(i).catIdx ==1 || tracks(i).catIdx ==2)
            AP2pos_idx=[AP2pos_idx idx(i)];
            tmp_pos_index=[tmp_pos_index i];
        end
    end
    % save labels for cmeviewer2_KM
    labels(1,AP2neg_idx)=2;
    labels(1,AP2neg_idx_long)=3;
    labels(1,AP2pos_idx)=5;
    temlabels{nd} = labels(Index);
    tracks_neg = tracks(tmp_neg_index);
    tracks_neg_long = tracks(tmp_neg_longindex);
    tracks_pos = tracks(tmp_pos_index);
    
    %% save AP2neg/AP2pos tracks information
    %AP2-neg
    if  length(tracks_neg)>=1
        tmp_tmp = getA_Lft_tdp(tracks_neg);
        [PLC_AP2neg_summary(nd).maxA,PLC_AP2neg_summary(nd).Lft,PLC_AP2neg_summary(nd).X,PLC_AP2neg_summary(nd).Y,PLC_AP2neg_summary(nd).Totaldisplacement,PLC_AP2neg_summary(nd).firstNmaxA,PLC_AP2neg_summary(nd).firstN_A] = ...
            deal(tmp_tmp.maxA,tmp_tmp.Lft,tmp_tmp.X,tmp_tmp.Y,tmp_tmp.Totaldisplacement,tmp_tmp.firstNframes_maxA,tmp_tmp.firstNframes);
        PLC_AP2neg_summary(nd).source = data(nd).source;
        PLC_AP2neg_summary(nd).number = length(tmp_tmp.maxA);
        PLC_AP2neg_summary(nd).cellArea = cellArea(nd);
        PLC_AP2neg_summary(nd).Coff_matrix = length(tmp_tmp.maxA)./(cellArea(nd).*data(nd).framerate.*data(nd).movieLength./60);
        PLC_AP2neg_summary(nd).mean_maxA = mean(tmp_tmp.maxA);
        PLC_AP2neg_summary(nd).mean_Lifetime = mean(tmp_tmp.Lft);
        PLC_AP2neg_summary(nd).mean_Totaldisplacement = mean(tmp_tmp.Totaldisplacement);
    end
    %AP2-neg-long
    if  length(tracks_neg_long)>=1
        tmp_tmp = getA_Lft_tdp(tracks_neg_long);
        [PLC_AP2negLong_summary(nd).maxA,PLC_AP2negLong_summary(nd).Lft,PLC_AP2negLong_summary(nd).X,PLC_AP2negLong_summary(nd).Y,PLC_AP2negLong_summary(nd).Totaldisplacement,PLC_AP2negLong_summary(nd).firstNmaxA,PLC_AP2negLong_summary(nd).firstN_A] = ...
            deal(tmp_tmp.maxA,tmp_tmp.Lft,tmp_tmp.X,tmp_tmp.Y,tmp_tmp.Totaldisplacement,tmp_tmp.firstNframes_maxA,tmp_tmp.firstNframes);
        PLC_AP2negLong_summary(nd).source = data(nd).source;
        PLC_AP2negLong_summary(nd).number = length(tmp_tmp.maxA);
        PLC_AP2negLong_summary(nd).cellArea = cellArea(nd);
        PLC_AP2negLong_summary(nd).Coff_matrix = length(tmp_tmp.maxA)./(cellArea(nd).*data(nd).framerate.*data(nd).movieLength./60);
        PLC_AP2negLong_summary(nd).mean_maxA = mean(tmp_tmp.maxA);
        PLC_AP2negLong_summary(nd).mean_Lifetime = mean(tmp_tmp.Lft);
        PLC_AP2negLong_summary(nd).mean_Totaldisplacement = mean(tmp_tmp.Totaldisplacement);
    end
    
    %AP2-pos
    if  length(tracks_pos)>=1
        tmp_tmp = getA_Lft_tdp(tracks_pos);
        [PLC_AP2pos_summary(nd).maxA,PLC_AP2pos_summary(nd).Lft,PLC_AP2pos_summary(nd).X,PLC_AP2pos_summary(nd).Y,PLC_AP2pos_summary(nd).Totaldisplacement,PLC_AP2pos_summary(nd).firstNmaxA,PLC_AP2pos_summary(nd).firstN_A] = ...
            deal(tmp_tmp.maxA,tmp_tmp.Lft,tmp_tmp.X,tmp_tmp.Y,tmp_tmp.Totaldisplacement,tmp_tmp.firstNframes_maxA,tmp_tmp.firstNframes);
        PLC_AP2pos_summary(nd).source = data(nd).source;
        PLC_AP2pos_summary(nd).number = length(tmp_tmp.maxA);
        PLC_AP2pos_summary(nd).cellArea = cellArea(nd);
        PLC_AP2pos_summary(nd).Coff_matrix = length(tmp_tmp.maxA)./(cellArea(nd).*data(nd).framerate.*data(nd).movieLength./60);
        PLC_AP2pos_summary(nd).mean_maxA = mean(tmp_tmp.maxA);
        PLC_AP2pos_summary(nd).mean_Lifetime = mean(tmp_tmp.Lft);
        PLC_AP2pos_summary(nd).mean_Totaldisplacement = mean(tmp_tmp.Totaldisplacement);
    end
end
PLC_summary.AP2_neg = PLC_AP2neg_summary;
PLC_summary.AP2_neg_long = PLC_AP2negLong_summary;
PLC_summary.AP2_pos = PLC_AP2pos_summary;

for iii = 1:length(data)
    if ~exist([data(iii).source 'Analysis'],'dir')
        mkdir([data(iii).source 'Analysis']);
    end
    labels = temlabels{iii};
    save([data(iii).source 'Analysis\labels.mat'], 'labels');
end
fprintf('filter done! \n');
filesepIndex = strfind(data(1).source,filesep);
rootpath = data(1).source(1:filesepIndex(end-3));
info = fopen([rootpath 'parameters.txt'],'a');
fprintf(info,'Filter Parameters\n');
fprintf(info,'Function plc_filterAP2neg_2ch\n');
fprintf(info,['date: ' datestr(now,31) '\n']);
fprintf(info,['Mask is [' num2str(ip.Results.mask) ']\n']);
fprintf(info,['A_SNR is [' num2str(asnr(1)) ' ' num2str(asnr(2)) ']\n']);
fclose(info);
save([rootpath 'PLC_infosummary.mat'],'PLC_summary');
end