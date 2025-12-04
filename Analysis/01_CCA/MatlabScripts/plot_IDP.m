function figH = plot_IDP(row_index)

load('/Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/IDP_MH_lookups.mat','IDP_lookup');
png_dir100 = '/Applications/workbench/UKBiobank_BrainImaging_GroupMeanTemplates/rfMRI_ICA_d100.sum/';
png_dir25 = '/Applications/workbench/UKBiobank_BrainImaging_GroupMeanTemplates/rfMRI_ICA_d25.sum/';

if strncmp(IDP_lookup.IDP_type(row_index),'AMP',3)
    node = IDP_lookup.Original_node{row_index};
    if strcmp(IDP_lookup.IDP_type(row_index),'AMP100')
        pic = imread([png_dir100 sprintf('%04d.png',node-1)]);
    elseif strcmp(IDP_lookup.IDP_type(row_index),'AMP25')
        pic = imread([png_dir25 sprintf('%04d.png',node-1)]);
    end
    figure; figH = imagesc(pic); axis off; axis equal;
else
    if strcmp(IDP_lookup.IDP_type(row_index),'FNET100') || strcmp(IDP_lookup.IDP_type(row_index),'PNET100')
        good = cell2mat(IDP_lookup.Original_node(1:55));
        node1 = good(IDP_lookup.Node1{row_index});
        pic1 = imread([png_dir100 sprintf('%04d.png',node1-1)]);
        node2 = good(IDP_lookup.Node2{row_index});
        pic2 = imread([png_dir100 sprintf('%04d.png',node2-1)]);
    elseif strcmp(IDP_lookup.IDP_type(row_index),'FNET25') || strcmp(IDP_lookup.IDP_type(row_index),'PNET25')
        good = cell2mat(IDP_lookup.Original_node(3026:3046));
        node1 = good(IDP_lookup.Node1{row_index});
        pic1 = imread([png_dir25 sprintf('%04d.png',node1-1)]);
        node2 = good(IDP_lookup.Node2{row_index});
        pic2 = imread([png_dir25 sprintf('%04d.png',node2-1)]);
    end
    figure
    subplot(1,2,1); figH.node1 = imagesc(pic1); axis off; axis equal;
    subplot(1,2,2); figH.node2 = imagesc(pic2); axis off; axis equal;
end
    

