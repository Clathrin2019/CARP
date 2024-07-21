function [] = filterCLTA(data, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', @isstruct);
ip.addParameter('mask',[],@isnumeric);
ip.addParameter('lifetime_cut',5 ,@isnumeric);
ip.addParameter('SNR',[1.05 1.05],@isnumeric);
ip.addParameter('A_SNR',[1.1 1.3],@isnumeric);
ip.addParameter('GAP_ratio',0.2,@isnumeric);
ip.addParameter('MaxIntensityThreshold', 1100,@isnumeric);
ip.parse(varargin{:});
snr = ip.Results.SNR;
lifetime_cut = ip.Results.lifetime_cut;
asnr = ip.Results.A_SNR;
gapratio = ip.Results.GAP_ratio;
intensity_cut = ip.Results.MaxIntensityThreshold;
CLTA_AP2neg_summary = struct;
CLTA_AP2pos_summary = struct;
cellArea = [];
parfor nd = 1:length(data)
    fprintf('filter start %d movies\n',nd);
    temtracks = load([data(nd).source 'Tracking\ProcessedTracks.mat']);
    tracks = temtracks.tracks;
    [mask,cellArea(nd)] = load_cellmask(data(nd),ip.Results.mask);   
    % sort lifetime for cmedataviewer2_KM
    [~,Index] = sort([tracks.lifetime_s], 'descend');
    labels = zeros(1,length(tracks));
    % prefilter
    AP2pos_idx=[];
    AP2neg_idx=[];
    [tmp_neg_index,tmp_pos_index]=deal([]); 
    for i = 1:length(tracks)
        Signalmch = tracks(i).A(1,:)+tracks(i).c(1,:);
        Signalmchlow = tracks(i).A(1,:)+tracks(i).c(1,:)-tracks(i).A_pstd(1,:);
        Noisehigh = tracks(i).c(1,:)+1.96*tracks(i).sigma_r(1,:);
        Noisemch = tracks(i).c(1,:);
        maxSignalmch = max(tracks(i).A(1,:)+tracks(i).c(1,:));
        tracks_meanSNR = mean(Signalmch)/mean(Noisemch);
        tracks_maxSNR = maxSignalmch/mean(Noisemch);
        x_position = round(nanmean(tracks(i).x(1,:)));
        y_position = round(nanmean(tracks(i).y(1,:)));
        gapnumber = length(find(tracks(i).gapVect==1));
        maskfilter = mask(y_position, x_position) == 1;
        catfilter = (tracks(i).catIdx==1);
        lifetimefilter = tracks(i).lifetime_s>=lifetime_cut;
        snrfilter = tracks_meanSNR>=snr(1) && tracks_maxSNR>=snr(2);
        a_defined_filter = mean(Signalmchlow)/mean(Noisehigh) >asnr(1);
        max_firstfilter = maxSignalmch/(tracks(i).A(1,1)+tracks(i).c(1,1)) > asnr(2);
        gapfilter = gapnumber/length(tracks(i).t) < gapratio;
        slavechfilter = (tracks(i).significantSlave(2)==1 && tracks(i).significantSlave(3)==1);
        slavechfilter2 = (tracks(i).significantSlave(2)==1 && tracks(i).significantSlave(3)==0);
        maxintensity_cut = (max(tracks(i).A(1,:)) > intensity_cut);
        if  maskfilter && lifetimefilter && catfilter && snrfilter && gapfilter && slavechfilter && maxintensity_cut
            AP2pos_idx = [AP2pos_idx i];
            tmp_pos_index = [tmp_pos_index i];
        end
        if  maskfilter && lifetimefilter && catfilter && snrfilter && gapfilter && slavechfilter2 
            AP2neg_idx = [AP2neg_idx i];
            tmp_neg_index = [tmp_neg_index i];
        end
    end
    % save labels for cmeviewer2_KM
    labels(1,AP2neg_idx)=2;
    labels(1,AP2pos_idx)=5;
    
    temlabels{nd} = labels(Index);
    tracks_neg = tracks(tmp_neg_index);
    tracks_pos = tracks(tmp_pos_index);
     %% save AP2neg/AP2pos tracks information
    %AP2-neg
    tmp_tmp = getA_Lft_tdp(tracks_neg);
    [CLTA_AP2neg_summary(nd).maxA,CLTA_AP2neg_summary(nd).Lft,CLTA_AP2neg_summary(nd).X,CLTA_AP2neg_summary(nd).Y,CLTA_AP2neg_summary(nd).Totaldisplacement,CLTA_AP2neg_summary(nd).firstNmaxA,CLTA_AP2neg_summary(nd).firstN_A] = ...
        deal(tmp_tmp.maxA,tmp_tmp.Lft,tmp_tmp.X,tmp_tmp.Y,tmp_tmp.Totaldisplacement,tmp_tmp.firstNframes_maxA,tmp_tmp.firstNframes);
    CLTA_AP2neg_summary(nd).source = data(nd).source;
    CLTA_AP2neg_summary(nd).number = length(tmp_tmp.maxA);
    CLTA_AP2neg_summary(nd).cellArea = cellArea(nd);
    CLTA_AP2neg_summary(nd).Coff_matrix = length(tmp_tmp.maxA)./(cellArea(nd).*data(nd).framerate.*data(nd).movieLength./60);
    CLTA_AP2neg_summary(nd).mean_maxA = mean(tmp_tmp.maxA);
    CLTA_AP2neg_summary(nd).mean_Lifetime = mean(tmp_tmp.Lft);
    CLTA_AP2neg_summary(nd).mean_Totaldisplacement = mean(tmp_tmp.Totaldisplacement);    
    %AP2-pos
    tmp_tmp = getA_Lft_tdp(tracks_pos);
    [CLTA_AP2pos_summary(nd).maxA,CLTA_AP2pos_summary(nd).Lft,CLTA_AP2pos_summary(nd).X,CLTA_AP2pos_summary(nd).Y,CLTA_AP2pos_summary(nd).Totaldisplacement,CLTA_AP2pos_summary(nd).firstNmaxA,CLTA_AP2pos_summary(nd).firstN_A] = ...
        deal(tmp_tmp.maxA,tmp_tmp.Lft,tmp_tmp.X,tmp_tmp.Y,tmp_tmp.Totaldisplacement,tmp_tmp.firstNframes_maxA,tmp_tmp.firstNframes);
    CLTA_AP2pos_summary(nd).source = data(nd).source;
    CLTA_AP2pos_summary(nd).number = length(tmp_tmp.maxA);
    CLTA_AP2pos_summary(nd).cellArea = cellArea(nd);
    CLTA_AP2pos_summary(nd).Coff_matrix = length(tmp_tmp.maxA)./(cellArea(nd).*data(nd).framerate.*data(nd).movieLength./60);
    CLTA_AP2pos_summary(nd).mean_maxA = mean(tmp_tmp.maxA);
    CLTA_AP2pos_summary(nd).mean_Lifetime = mean(tmp_tmp.Lft);
    CLTA_AP2pos_summary(nd).mean_Totaldisplacement = mean(tmp_tmp.Totaldisplacement);
    
end
CLTA_summary.AP2_neg = CLTA_AP2neg_summary;
CLTA_summary.AP2_pos = CLTA_AP2pos_summary;
%%
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
fprintf(info,'Function filterCLTA\n');
fprintf(info,['date: ' datestr(now,31) '\n']);
fprintf(info,['Mask is [' num2str(ip.Results.mask) ']\n']);
fprintf(info,['SNR is [' num2str(asnr(1)) ' ' num2str(asnr(2)) ']\n']);
fclose(info);
save([rootpath 'CLTA_infosummary.mat'],'CLTA_summary');
end

