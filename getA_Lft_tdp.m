function [statisticalData] = getA_Lft_tdp(tracks)
for i=1:length(tracks)
    if ~isempty(tracks)
        statisticalData.maxA(i,1) = nanmax(tracks(i).A(1,:));
        if length(tracks(i).A(1,:))>=15
            statisticalData.firstNframes_maxA(i,1) = nanmax(tracks(i).A(1,1:15));
            statisticalData.firstNframes{i} = tracks(i).A(1,1:15);
        else
            statisticalData.firstNframes_maxA(i,1) = NaN;
            statisticalData.firstNframes{i} = [];
        end
        statisticalData.Lft(i,1) = tracks(i).lifetime_s;
        statisticalData.X(i,1) = mean(tracks(i).x(1,:));
        statisticalData.Y(i,1) = mean(tracks(i).y(1,:));
        
        if ~isempty(tracks(i).MotionAnalysis)
            statisticalData.Totaldisplacement(i,1) = tracks(i).MotionAnalysis.totalDisplacement;
        end
    else
        statisticalData.maxA(i,1)=[];
        statisticalData.Lft(i,1)=[];
        statisticalData.Totaldisplacement(i,1)=[];
        statisticalData.X(i,1) = [];
        statisticalData.Y(i,1) = [];
    end
end
end

