% 2025 09 21
beep off

% This script will be run on the CARC to:
%   Read in the large data file,
%   Extract the data per sample,
%   Cluster the individual mets, and
%   Write each met to an individual file.
%   Visualize each met, in the sample, with clustering.

% Output will be the files that contain the identified clusters, with the
% x, y, cell state and phenotypes.

% TK 10.14.25 -- changing output to include all markers

saveImages = 1;
writeToFile = 1;

m0 = hex2rgb('#d6d6d6');  % m0 grey
m1 = hex2rgb('#88e0b0'); % anti tumor Tam - light green
m2 = hex2rgb('#2fa130'); % Pro tumor tam Dark green
cancer = hex2rgb('#2278b3'); % Cancer - blue
cd4_th = hex2rgb('#fc9a97'); % CD4 Light Red
treg = hex2rgb('#e11d18'); % Treg - Red
cd8 = hex2rgb('#febe6f'); % CD8 - Orange
nk = hex2rgb('#c7b3d3');  % NK - Lavender
mdsc = hex2rgb('#f5e205');  % MDSC - Yellow
     
apriori_sample_name = ["121223 P9msP61-2 #01 NT2.5 LM Lungs 10d-1";
"121223 P9msP61-2 #02 NT2.5 LM Lungs 10d-2";
"121223 P9msP61-2 #03 NT2.5 LM Lungs 10d-3";
"121223 P9msP61-2 #04 NT2.5 LM Lungs 10d-4";
"121223 P9msP61-2 #05 NT2.5 LM Lungs 10d-5";
"121223 P9msP61-2 #06 NT2.5 LM Lungs 16d-1";
"121223 P9msP61-2 #08 NT2.5 LM Lungs 16d-3";
"121223 P9msP61-2 #10 NT2.5 LM Lungs 22d-1";
"121223 P9msP61-2 #12 NT2.5 LM Lungs 22d-3";
"121223 P9msP61-2 #13 NT2.5 LM Lungs 22d-4";
"121223 P9msP61-2 #14 NT2.5 LM Lungs 22d-5";
"121223 P9msP61-2 #15 NT2.5 LM d122 UN TX 1";
"121223 P9msP61-2 #16 NT2.5 LM d122 V1";
"121223 P9msP61-2 #17 NT2.5 LM day33 V8";
"121223 P9msP61-2 #19 7 days Untx-1";
"121223 P9msP61-2 #20 7 days Untx-2";
"121223 P9msP61-2 #21 7 days Untx-3";
"121223 P9msP61-2 #22 7 days Untx-4";
"121223 P9msP61-2 #23 NeuN NT2.5 LM 5 wks #4";
"121223 P9msP61-2 #24 NeuN NT2.5 LM 5 wks #5";
"121223 P9msP61-2 #25 NeuN NT2.5 LM 5 wks #1";
"121223 P9msP61-2 #27 NeuN NT2.5 LM 5 wks #3";
"121223 P9msP61-2 #28 NeuN NT2.5 LM 4 wks #1";
"121223 P9msP61-2 #29 NeuN NT2.5 LM 4 wks #2";
"121223 P9msP61-2 #30 NeuN NT2.5 LM 4 wks #3";
"121223 P9msP61-2 #31 NeuN NT2.5 LM 4 wks #4";
"121223 P9msP61-2 #32 NeuN NT2.5 LM 4 wks #5"];
 
 apriori_num_mets_per_sample =[2, 2, 3, 2, 1, 2, 1, 1, 1, 3, 3, 11, 2, 4, 7, ...
     1, 6, 4, 5, 4, 4, 4, 16, 15, 3,12, 18];
 
% pathways
imaging_data_pathway = "./data/WholeTissueConsolidated_data.csv";
output_pathway = "./output/";
sample_files_pathway = strcat(output_pathway,"samples/");
cluster_output_pathway = strcat(output_pathway,"metastases/");

wholeTissue = readtable(imaging_data_pathway); % import data
unique_ROIs = unique(wholeTissue.SampleName); % find unique ROI labels
sample_name = unique(extractBefore(unique_ROIs,"_[")); % find unique sample names

% This table has 9 columns: the phenotypes from the mIHC and a cell type
% label (in text)
cell_definition_table = readtable('./data/definitions.xlsx','TextType','string'); 
cell_state_values = [0:6 8 10 -10 -10];
cell_state_labels = ["M0","M1","M2","cancer","CD4","treg","CD8","NK","MDSC","immune","unknown"];

col_names = ["Name","Unclassified","M1","M2","Cancer","CD4","Treg","CD8","NK","MDSC","Total"];
col_type = {'string','double','double','double','double','double','double','double','double','double','double'};

table_idx = 0;

%% Build the cell definition look up table

% Clean up table
for c = 1:9
    cell_definition_table.(c) = erase(string(cell_definition_table.(c)), '''');
end

% Extract marker phenotypes and labels
phenotypes       = cell_definition_table{:, 1:8}; % 8 marker columns (string now)
phenotype_labels = cell_definition_table{:, 9};   % cell type names (string)

% Convert labels to cell_state values
phenotype_code = zeros(height(cell_definition_table), 1);
for i = 1:height(cell_definition_table)
    match_idx = find(strcmp(phenotype_labels(i), cell_state_labels));
    if ~isempty(match_idx)
        phenotype_code(i) = cell_state_values(match_idx);
    else
        phenotype_code(i) = -10; % handle missing/unknown
    end
end

% Build keys for marker patterns
keys = strings(height(cell_definition_table), 1); % preallocate string array
for i = 1:height(cell_definition_table)
    row = string(phenotypes(i, :));     % ensure row is string array
    keys(i) = strjoin(row, '_');        % scalar string key
end

% Make map from pattern to cell state value
map = containers.Map(keys, phenotype_code);



%% Process each sample, write data to ROI specific files.
for i = 1:numel(sample_name)

    % Find data associated with that sample
    logicalArray = contains(wholeTissue.SampleName,sample_name{i});
    sample_data = wholeTissue(logicalArray, 1:246); % This eliminates too much data; keep all of it
    cell_state = -10*ones(height(sample_data),1);

    % X -> 6, Y -> 7, phenos -> 239:246
    for j = 1:height(sample_data)
        row_markers = string(sample_data{j, 239:end});
        
        key = strjoin(row_markers,'_');

        if isKey(map,key)
            cell_state(j) = map(key);
        else 
            cell_state(j) = -10;
        end
    end

    sample_data = addvars(sample_data, cell_state,'After','CellYPosition','NewVariableNames','cell_state'); % adds sample data after CellYPosition (pos 8)

    filePathway = strcat(sample_files_pathway,'dataFiles/',sample_name{i},".csv");
    if writeToFile
        writetable(sample_data,filePathway);
    end

    % Clustering:
    % find the corresponding number of mets in the apriori array,
    idx = find(contains(apriori_sample_name,sample_name{i}));
    if ~isempty(idx)
        
        T = clusterdata(table2array(sample_data(:, 6:7)),'linkage','single','maxclust',apriori_num_mets_per_sample(idx)); % X -> 6, Y -> 7 in table!
    
        h = figure('units','normalized','outerposition',[0 0 0.5 0.7]);
        cmap = colormap(winter(apriori_num_mets_per_sample(idx)));
        t= tiledlayout('flow');
        ax(1)=nexttile;
        gscatter(table2array(sample_data(:,6)), table2array(sample_data(:,7)), T,cmap); % X -> 6, Y -> 7
        lg = legend('Location','northoutside','Orientation','horizontal');
        lg.Box = "off";
        yl = ylim;
        xl = xlim;
        axis square
        box on
        fontsize(16,"pixels");
    
        for j = 1:apriori_num_mets_per_sample(idx)
            table_idx = table_idx + 1;
            ax(j+1)=nexttile;
            clusterID = find(T==j);
            clusterData = sample_data(clusterID,:);
    
            xcenter = (min(table2array(clusterData(:,6))) + max(table2array(clusterData(:,6))))/2;
            ycenter = (min(table2array(clusterData(:,7))) + max(table2array(clusterData(:,7))))/2;
            shiftedData = [table2array(clusterData(:,6))-xcenter table2array(clusterData(:,7))-ycenter];
            s=scatter(shiftedData(:,1),shiftedData(:,2),20,'.','MarkerEdgeColor',cmap(j,:));
            s.SizeData = 10;
            axis square
            box on
            title(strcat("Cluster ", num2str(j)))
    
            % Save the data per met; shifted X & Y go first
            if writeToFile
                writetable([array2table(shiftedData) sample_data(clusterID,1:end)],strcat(cluster_output_pathway,"/dataFiles/",sample_name{i},"_met_",num2str(j),".csv"));
            end
        end
    
        set(ax,'XTickLabels',[])
        set(ax,'YTickLabels',[])
        title(t,sample_name{i},'Interpreter', 'none')

        if saveImages
            exportgraphics(h,strcat(sample_files_pathway,"Images/",sample_name{i},"_cluster.jpg"));
            savefig(h,strcat(sample_files_pathway,"Images/",sample_name{i},"_cluster.fig"));
        end
    end
end

