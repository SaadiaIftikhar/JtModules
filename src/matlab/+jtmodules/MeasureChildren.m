classdef  MeasureChildren < handle

properties (Constant)

    VERSION = '0.0.3'
end

methods (Static)    

function [measuremnets, figure] =  main(matImpSpotFeatures, matParentCount, ...
         matParentObject, plot)


% Help for the Measure Texture module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Pools measurements of children to create measurments for each parent. The
% measurments include mean/median and var as well as higher central
% moments. Values with NaNs are excluded from the analysis.
% *************************************************************************
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
%%% DATA ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%

% drawnow

numParents = max([max(matParentObject(:)), matParentCount]);

% initialize output
PCNaNMean = nan(numParents,size(matImpSpotFeatures,2));
PCNaNMedian = nan(numParents,size(matImpSpotFeatures,2));
PCNaNVar = nan(numParents,size(matImpSpotFeatures,2));
PCNaNMoment3 = nan(numParents,size(matImpSpotFeatures,2));
PCNaNMoment4 = nan(numParents,size(matImpSpotFeatures,2));
PCNaNMoment5 = nan(numParents,size(matImpSpotFeatures,2));
PCNaNMoment6 = nan(numParents,size(matImpSpotFeatures,2));
PCNaNStd = nan(numParents,size(matImpSpotFeatures,2));


if any(matParentObject>0) && (numParents>0);
    
    % sort measurments of children according to parents. this will speed up
    % the later calculation
    [smatParentObject sIX] = sort(matParentObject); % sort children according to their parents
    
    smatImpSpotFeatures = matImpSpotFeatures(sIX,:); % sort data of children so that data of siblings is next to each other. This allows to use a fast strategy for filtering
    
    [usmatParentObject bFirst] = unique(smatParentObject,'first'); % get ID of first child
    [~, bLast] = unique(smatParentObject,'last'); % last child
    
    for j=1:numParents % for each parent
        f = usmatParentObject == j; % identify children
        if any(f)
            CurrData = smatImpSpotFeatures(bFirst(f):bLast(f),:); % import data of all children
            
            PCNaNMean(j,:) =        nanmean(CurrData,1);
            PCNaNMedian(j,:) =      nanmedian(CurrData,1);
            PCNaNVar(j,:) =         nanvar(CurrData,[],1);
            PCNaNStd(j,:) =         nanstd(CurrData,[],1);

            
            for k=1:size(smatImpSpotFeatures,2) % loop through each feature. note that it is possible that only some features have nans.
                
                ff = ~(isnan(CurrData(:,k)));
                CurrColumnData = CurrData(ff,k);
                
                PCNaNMoment3(j,k) =     moment(CurrColumnData,3,1);
                PCNaNMoment4(j,k) =     moment(CurrColumnData,4,1);
                PCNaNMoment5(j,k) =     moment(CurrColumnData,5,1);
                PCNaNMoment6(j,k) =     moment(CurrColumnData,6,1);

            end
        end
    end
end

measuremnets(1) = PCNaNMean; 
measuremnets(2) = PCNaNMedian; 
measuremnets(3) = PCNaNMoment3;
measuremnets(4) = PCNaNMoment4;
measuremnets(5) = PCNaNMoment5;
measuremnets(6) = PCNaNMoment6; 
measuremnets(7) = PCNaNStd;


%%%%%%%%%%%%%%%%%%%%
%%% SAVE RESULTS %%%
%%%%%%%%%%%%%%%%%%%%


% note that for each parent derived measuremnt is saved in new file. 

figure = '';

end
end
end
