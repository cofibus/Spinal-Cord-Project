classdef SpinalCordGenetics 
    methods (Static) 
        function Genetics = GetMouseGenetics(Sym) 
            load(['03_SpinalCordGeneticsFiles' filesep 'matCountSpatialScaled.mat'])
            load(['03_SpinalCordGeneticsFiles' filesep 'matCountSingleCell.mat'])
            load(['03_SpinalCordGeneticsFiles' filesep 'matRCTD.mat']) 
            load(['01_SpinalGeometryFiles' filesep 'matGeometry.mat']);
            load(['03_SpinalCordGeneticsFiles' filesep 'matCellInfo.mat'])
            load(['03_SpinalCordGeneticsFiles' filesep 'matNeuronChat.mat'])
            Genetics.Counts = CountSpatial;
            Genetics.Genes = GenesSpatial;
            Genetics.Weights = weights;
            Genetics.Types = RCTDTypes;    
            Genetics.Coords = RCTDCoord;
            Genetics.Census = CellTypesCensus;
            Genetics.CountsSC = CountSingleCell;
            Genetics.TypesSC = Types;
            Genetics.GenesSC = Genes;
            Genetics.LineageSC = Lineage;

            Genetics.NCTypes = NCTypes;
            clus = categorical(erase(string(Genetics.Census.Cluster),'-'));
            for ii = 1:length(Genetics.NCTypes)
               Genetics.NCTypes(ii) = Genetics.Census.Type(clus == Genetics.NCTypes(ii));
            end
            Genetics.NCProba = NCProba';
            %Genetics.NCProba(:) = Genetics.NCProba(:)./max(Genetics.NCProba(:));
            Genetics.NCProba(Genetics.NCProba~=0 & ~isnan(Genetics.NCProba)) = normcdf(zscore(Genetics.NCProba(Genetics.NCProba~=0  & ~isnan(Genetics.NCProba))));

            unq = unique(Genetics.NCTypes,'stable');
            ncp = [];
            for ii = 1:length(unq)
                whr = ismember(Genetics.NCTypes,unq(ii));
                ncp(ii,:) = mean(Genetics.NCProba(whr,:),1,'omitnan');
            end
            ncp2 =[];
            for ii = 1:length(unq)
                whr = ismember(Genetics.NCTypes,unq(ii));
                ncp2(:,ii) = mean(ncp(:,whr),2,'omitnan');
            end
            Genetics.NCProba = ncp2;
            Genetics.NCProba(isnan(Genetics.NCProba)) = 1;
            Genetics.NCTypes = unq;
        end
    end
end