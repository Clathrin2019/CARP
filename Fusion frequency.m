function [] = Fusion frequency(varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('AP2ch',2, @isnumeric); % 2-561nm 3-637/647nm
ip.addParameter('usefulChNum',2, @isnumeric); % 2 channels or 3 channels
ip.addParameter('Parameters', [1.49 200 16e-6], @(x) numel(x)==3);
ip.addParameter('filter', true, @islogical);
ip.parse(varargin{:});
parameters = ip.Results.Parameters;
AP2ch = ip.Results.AP2ch;
usefulChNum = ip.Results.usefulChNum;
data_path = uigetdir('Z:\Data_He Lab');
if usefulChNum == 3
    data = loadConditionData(data_path,{'EGFP', 'TagRFP', 'JF646'},'markers',{'gfp','rfp','Cy5'},'Parameters', parameters);
    if ip.Results.filter
        filterAP2neg_3ch(data,'SNR',[1.05 1.6],'mask',[]);
    end
else
    if AP2ch==2
        data = loadConditionData(data_path,{'EGFP', 'TagRFP'},'markers',{'gfp','rfp'},'Parameters', parameters);
        if ip.Results.filter
            filterAP2neg_2ch(data,'SNR',[1.05 1.6],'mask',[]);
        end
    else
        data = loadConditionData(data_path,{'EGFP', 'JF646'},'markers',{'gfp','Cy5'},'Parameters', parameters);
        if ip.Results.filter
            filterAP2neg_2ch(data,'SNR',[1.05 1.6],'mask',[],'AP2ch',3);
        end
    end
end
tmp_path = dir(data_path);
subfolders_path = [];
k = 1;
%%
for i =1:length(tmp_path)
    if tmp_path(i).isdir && ~isequal(tmp_path(i).name, '.') && ~isequal(tmp_path(i).name, '..')
        subfolders_path{k,1} = [tmp_path(i).folder filesep tmp_path(i).name];
        protein_name{k,1} = tmp_path(i).name;
        k = k + 1;
    end
end
%%

parfor i =1:length(subfolders_path)
    if AP2ch==2
        tmpdata = loadConditionData(subfolders_path{i},{'EGFP', 'TagRFP'},'markers',{'gfp','rfp'},'Parameters', parameters);
    else
        tmpdata = loadConditionData(subfolders_path{i},{'EGFP', 'JF646'},'markers',{'gfp','Cy5'},'Parameters', parameters);
    end
    AP2_negcount = [];
    cellArea = [];
    tmpAP2_negLft = [];
    tmpAP2_A = [];
    for j =1:length(tmpdata)
        labels = load([tmpdata(j).source '\Analysis\labels.mat']);
        tracks = load([tmpdata(j).source '\Tracking\ProcessedTracks.mat']);
        labels = labels.labels;
        tracks = tracks.tracks;
        [~, sort_index] = sort([tracks.lifetime_s], 'descend');
        label2=labels==2;
        label2 = sort_index(label2);
        tracks2 = tracks(label2);
        [~,cellArea(j)] = load_cellmask(tmpdata(j),[]);
        AP2_negcount(j,1) = length(find(labels==2));
        if AP2_negcount(j,1) ~=0
            tmpAP2_negLft(j,1) = mean([tracks2.lifetime_s]);
            for kk = 1:length(tracks2)
                tracks2(kk).maxA = max(tracks2(kk).A(1,:));
            end
            tmpAP2_A(j,1) = mean([tracks2.maxA]);
        else
            tmpAP2_negLft(j,1) = NaN;
            tmpAP2_A(j,1) = NaN;
        end
    end
    tmp_AP2_negCoff_matrix{i} = AP2_negcount./(cellArea.*tmpdata(1).framerate.*tmpdata(1).movieLength./60)';
    tmp_AP2_negLft{i} = tmpAP2_negLft;
    cell_Area{i} = cellArea;
    tmp_AP2_neg_A{i} = tmpAP2_A;
end

num_data = cellfun(@length,tmp_AP2_negLft);
infoStruct = [];
for m = 1:length(num_data) 
    for n = 1:num_data(m)
        infoStruct(end+1).name = protein_name{m};
        infoStruct(end).AP2neg_Lft = tmp_AP2_negLft{m}(n);
        infoStruct(end).AP2neg_coeff = tmp_AP2_negCoff_matrix{m}(n);
        infoStruct(end).AP2neg_A = tmp_AP2_neg_A{m}(n);
        infoStruct(end).cellArea = cell_Area{m}(n);
    end
end
infoTable = struct2table(infoStruct);
writetable(infoTable,[data_path '\plc_summary.csv']);
end