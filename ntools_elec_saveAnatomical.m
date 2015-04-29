function [anatomical_text, EOI] = ntools_elec_saveAnatomical(subj,hemi,elec_bin,elec_text)

% paint the electrodes onto subject's pial surface and output the
% anatomical regions (in percentage) where each electrode locates 
%
% Usage: ntools_elec_saveAnatomical(subj,hemi,elec_text)
% default cortical parcellation: ?h.aparc.annot
%
% Input:
% subj: subject ID in SUBJECTS_DIR
% hemi: hemisphere
% elec_text: subject electrode location file in T1 space
% elec_bin: subject electrode nifti image
%
% Output:
% cortical_text: anatomical regions of each electrode in percentage
%
% created by Hugh Wang, 3/11/2015, Xiuyuan.Wang@nyumc.org
%

%% read in the text file and parse the G/S/D electrodes

fprintf('\n%s\n',subj);

fid = fopen(elec_text);
elec_all = textscan(fid,'%s %f %f %f %s');
elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];

% Separate Grid, Strip and Depth electrodes

if isempty(char(elec_all{5}(:)))
%     g = strncmpi('G',elec_cell(:,1),1);
    d = strncmpi('D',elec_cell(:,1),1);
else
%     g = strncmpi('G',elec_all{5},1);
    d = strncmpi('D',elec_all{5},1);
end

elec_gs = elec_cell(~d,:);
elec_depth = elec_cell(d,:);


%% process with G/S

PathName = fileparts(elec_text);
cfg = [];
cfg.outdir = [PathName '/labels'];
    
hippo_elec = cell(1);
entorhinal_elec = cell(1);
lateral_frontal_elec = cell(1);

if ~isempty(elec_gs)
    % load pial surface
    [surf] = fs_load_subj(subj,hemi,'pial');
    [surf] = fs_calc_triarea(surf);
    surf.coords = surf.vertices;
    
    PathName = fileparts(elec_text);
    cfg = [];
    cfg.subject = subj;
    cfg.surf= surf;
    cfg.elec_names = elec_gs(:,1);
    cfg.elec_coords = cell2mat(elec_gs(:,2:4));
    cfg.hemi = hemi;
    cfg.outdir = [PathName '/labels'];
    cfg.fsavg = 0;
    
    annotfile = ntools_elec_saveAnnot(cfg);
    
    %% loading ?h.aparc.annot and get the region percentage of each elec
    
    % read annotation
    [~, elec_label, elec_colortable] = fs_read_annotation(annotfile);
    
    [~, label, colortable] = fs_read_annotation([getenv('SUBJECTS_DIR'),subj,'/label/',hemi,'.aparc.annot']);
    
    
    %%   
    for i=1:length(elec_gs)
        note = [];
        elec_nbrs = elec_label==elec_colortable.table(i,5);
        
        % didn't find neighbours in the label, happens when the electrode share
        % the same location with others
        if sum(elec_nbrs)==0
            elec_gs(i,5) = {[]};
            continue;
        end
        
        total_area = surf.vertex_area(elec_nbrs);
        aparc_table_idx = label(elec_nbrs);
        
        % remove 0 index
        aparc_table_idx = aparc_table_idx(aparc_table_idx>0);
        if isempty(aparc_table_idx), continue; end;
        
        aparc_idx_uniq = unique(aparc_table_idx);
        for j=1:length(aparc_idx_uniq)
            aparc_region = colortable.struct_names{colortable.table(:,5)==aparc_idx_uniq(j)};
            aparc_area = sum(total_area(aparc_table_idx==aparc_idx_uniq(j)));
            aparc_area_ratio = aparc_area/sum(total_area)*100;
            
            note = [note, sprintf('%0.2f%% %s ',aparc_area_ratio,aparc_region)];
        end
        
        if sum(strcmp('parahippocampal',strsplit(note,' ')))
            hippo_elec(end+1) = elec_gs(i,1);
        end
        if sum(strcmp('entorhinal',strsplit(note,' ')))
            entorhinal_elec(end+1) = elec_gs(i,1);
        end
        if sum(strcmp('caudalmiddlefrontal',strsplit(note,' '))) || ...
                sum(strcmp('rostralmiddlefrontal',strsplit(note,' '))) || ...
                sum(strcmp('parsopercularis',strsplit(note,' '))) || ...
                sum(strcmp('parstriangularis',strsplit(note,' '))) || ...
                sum(strcmp('parsorbitalis',strsplit(note,' '))) || ...
                sum(strcmp('lateralorbitofrontal',strsplit(note,' '))) || ...
                sum(strcmp('frontalpole',strsplit(note,' ')))
            lateral_frontal_elec(end+1) = elec_gs(i,1);
        end
        
        elec_gs(i,5) = cellstr(note);
        clear elec_nbrs total_area aparc*
    end

end

%% process depth electrodes 
if ~isempty(elec_depth)
    hdr = ntools_elec_load_nifti(elec_bin);
    depth_row = find(d);

    % convert aseg
    aseg_mgz = [getenv('SUBJECTS_DIR'),subj,'/mri/aseg.mgz'];
    aseg_nii = [cfg.outdir,'/aseg.nii.gz'];
    if ~exist(aseg_nii,'file')
        [status,msg] = unix(sprintf('mri_convert --out_orientation LAS %s %s',aseg_mgz,aseg_nii));
        if status, disp(msg); return; end;
    end
    aseg = ntools_elec_load_nifti(aseg_nii);

    [seg_idx, seg_name] = xlsread('/space/mdeh1/5/halgdev/projects/nyuproj/loc/Milan/aparc_aseg_idx_name.xlsx');

    for k=1:length(depth_row)
        note = [];
        
        seg_num = aseg.vol(hdr.vol==depth_row(k));
        
        % continue if seg_num is empty
        if isempty(seg_num), elec_depth(k,5) = {[]}; continue; end;
        
        unique_seg_num = unique(seg_num);
        for m=1:length(unique_seg_num)
            seg_num_vox = sum(seg_num==unique_seg_num(m));
            note = [note, sprintf('%.2f%% %s ',100*seg_num_vox/27, seg_name{seg_idx==unique_seg_num(m)})];
        end
        
        if sum(strcmp('Left-Hippocampus',strsplit(note,' '))) || ...
           sum(strcmp('Right-Hippocampus',strsplit(note,' ')))
            hippo_elec(end+1) = elec_depth(k,1);
        end
        elec_depth(k,5) = cellstr(note);
    end    
end

%% save into text file

hippo_elec{1} = length(hippo_elec)-1;
entorhinal_elec{1} = length(entorhinal_elec)-1;
lateral_frontal_elec{1} = length(lateral_frontal_elec)-1;

anatomical_text = [PathName,'/',subj,'_T1_',hemi,'_AnatomicalRegions.txt'];
ntools_elec_savetxt(anatomical_text,[elec_gs;elec_depth]);

fid = fopen(anatomical_text,'a');
fprintf(fid,'\n');
fprintf(fid,'%% Total number of electrodes %d\n',length(elec_cell));

fprintf(fid,'%% Number of hippocampus electrodes %d\n',hippo_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(hippo_elec), fprintf(fid,'%s ',hippo_elec{ll}); end; fprintf(fid,'\n');

fprintf(fid,'%% Number of entorhinal cortex electrodes %d\n',entorhinal_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(entorhinal_elec), fprintf(fid,'%s ',entorhinal_elec{ll}); end; fprintf(fid,'\n');

fprintf(fid,'%% Number of lateral frontal electrodes %d\n',lateral_frontal_elec{1});
fprintf(fid,'%% ');
for ll = 2:length(lateral_frontal_elec), fprintf(fid,'%s ',lateral_frontal_elec{ll}); end; fprintf(fid,'\n');

fclose(fid);


EOI = [length(elec_cell),hippo_elec{1},entorhinal_elec{1},lateral_frontal_elec{1}];



