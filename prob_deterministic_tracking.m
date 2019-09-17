function prob_deterministic_tracking(dwi_folder,out_dir,nsim,atlas4connectome,t1_folder,preproc_func_handle2)
%function is in a weird state because I'm trying to make this more modular

% NOTES FOR RUNNING:
%{
%1
%d2n2s (repository) must be on path, clone from https://github.com/andrew-taylor-2/d2n2s
%and put on matlab path

%2
%I tried to cut out all the steps I did that didn't work. I ended up trying stuff on several different directories, so hopefully all my paths work here.
%Hopefully the script still flows.

%3 
%certain paths need to be changed. I think I put "%CHANGE" after each one.
%also at one point you'll have to give a mask you made yourself. also you
%have to tell dke ur number of directions in a variable at one point.

%4
%while it shouldn't be an issue in the body of the script, in the first two
%"d2n2s" commands below, you have to give d2n2s a folder with ONLY a .nii,
%a .bvec, a .bval, and optionally a .json. If the folder contains multiple
%files of these types, it won't work right. Also, the files need to match
%in number of elements and stuff. You can't have a .nii file with 61 images
%and a .bval or .bvec with only 60 lines.
%}



%% give anonymous functions to everyone
% is it really necessary tho

paren=@(x,varargin) x(varargin{:});
curly=@(x,varargin) x{varargin{:}};
append=@(x,str) fullfile(choose_output(@() fileparts(x),1),[choose_output(@() fileparts(x),2) str choose_output(@() fileparts(x),3)]); %changed slightly so that you don't need path
prepend=@(x,str) fullfile(choose_output(@() fileparts(x),1),[str choose_output(@() fileparts(x),2) choose_output(@() fileparts(x),3)]); %changed slightly so that you don't need path

do.paren=paren;
do.curly=curly;
do.appendd=append;
do.prepend=prepend;


%% sanitize inputs and give errors early if at all

if ~exist('atlas4connectome','var') || isempty(atlas4connectome) || ~exist(atlas4connectome,'file')
    warning('can''t find connectome file. Making tracts but not connectome')
end

%make sure out_dir makes sense. if it doesn't exist, it'll be made when
%images are written
assert(exist('out_dir','var'))
out_dir=char(out_dir);
[a,b,~]=fileparts(out_dir); %this does preclude using folders that contain a period...
out_dir=fullfile(a,b);

%% MAIN

[folder_in_which_to_write,working_dwi_name,dwi_first_b0]=denoise_and_create_sims(dwi_folder,out_dir,nsim,do)

[rt1,gmwmi_mask,mask_name,out_5tt]=t1_preproc(t1_folder,dwi_first_b0,out_dir)


[folder_in_which_to_write,working_dwi_name]=preproc_func_handle(folder_in_which_to_write,working_dwi_name,do) % this is kind of a stand in for what the pydesigner pipeline will eventually be
% if you wanted to use this as is, you could make a wrapper for py designer
% that just parses input and output file names and hands them to system.
% Also, here preproc_func_handle isn't an input, but it could be.

[]=track_each(folder_in_which_to_write,working_dwi_name,mask_name,gmwmi_mask)

end






function [folder_in_which_to_write,dwi_first_b0,working_dwi_name]=denoise_and_create_sims(dwi_folder,out_dir,nsim,do)

%% read DWIs
% FUTURE -- add support for multiple input folders

%---
%DWI
%---

%read
ff.no='bvecbvalimgjson';
dwi_init=d2n2s(dwi_folder,ff); %kind of a roundabout way of getting a filename but whatever

%-----
%NOISE
%-----

% get noise
[out1,out2]=mrtrix_denoise(dwi_init(1).fns.nii,'_dn') %blanks are shell and exe. no need for outs really
%out names
noise_out_name=fullfile(dwi_folder,'noise.nii');
denoised_dwi_name=do.append(dwi_init(1).fns.nii,'_dn');

%get the noise img
noise_img=spm_read_vols(spm_vol(noise_out_name)); % i can justify d2n2s if i'm at least wrapping a dir call but not if i'm just loading an image lawl

%constraint and reshape
noise_img(isnan(noise_img))=0;
noise_img=reshape(noise_img,[],1);

%------
%DWI_DN
%------
ffdn.pick=denoised_dwi_name;
dwi_dn=d2n2s(dwi_folder,ffdn);

%sizes
nvol=numel(dwi_dn);
dwi_size=size(dwi_dn(1).img);

%reshape
dwi_dn_img=cat(4,dwi_dn.img);
dwi_dn_img=reshape(dwi_dn_img,[],nvol);

%grab the first for t1 alignment
dwi_first_b0=dwi_dn(do.paren(find([dwi_dn.bval]<101),1)); %grabs the first b0


%% begin DWI sim and processing

%make a big ole matrix of random numbers
randomm=rand(prod(dwi_size),nvol,nsim,2);

%give a working dwi name just to make things a little more programmatic
working_dwi_name='dwi';

%% create and write realizations
tic

noisefunc=@(im_noise,rand_img1,rand_img2) sqrt(2.*im_noise.^2).*(erfinv(2.*rand_img1-1) + 1i.*erfinv(2.*rand_img2-1));
sim_noise_img=noisefunc(noise_img,randomm(:,:,:,1),randomm(:,:,:,2));

signal=abs(dwi_dn_img+sim_noise_img);%signal is now all blowed out in the three d: voxels by volumes by nsims
signal_img = reshape(signal,[dwi_size,nvol,nsim]);
signal_img = squeeze(num2cell(signal_img,[1 2 3])); %reassigning to same variable name to keep memory usage lower

folder_in_which_to_write=cell(nsim,1);
for u = 1:nsim %this could probably be parallelized -- but passing signal_imgc into the loop would be annoying overhead. Should just make a variable for each segment I think.
    
    fprintf('\n simulation andrew ver #%03d \n',u)

    [dwi_dn.img]=signal_img{:,u};
    
    f.gr=0; %this was necessary to put into designer, but you could use f=[] if you wanted
    folder_in_which_to_write{u}=fullfile(out_dir,['sim' num2str(u)]);
    d2n2s_write(dwi_dn,folder_in_which_to_write{u},working_dwi_name,f)    
    
end

reall_time=toc;
fprintf('Time to add noise and write realizations: %s',num2str(reall_time))

end



function [folder_in_which_to_write,working_dwi_name]=preproc_func_handle(folder_in_which_to_write,working_dwi_name,do)
%% I forgot to do this before, but ya gotta denoise it again
% I can just put it in the folder it was in before and use pick later

%in and out name for denoising
in=[working_dwi_name '.nii'];
out=do.append(in,'_dn');
[~,working_dwi_name,~]=fileparts(out);

%denoise and get the output names
dn_and_give_out_name=@(x) do.curly({system(['dwidenoise ' fullfile(x,in) ' ' fullfile(x,out)]) , fullfile(x,out)},2);
dwis_preprocd=cellfun(dn_and_give_out_name,folder_in_which_to_write,'un',0);

%just grab these and write to the out folder (the use of d2n2s is justified
%here IMO bc it wraps 4 dir calls and 4 writes)
for i=1:numel(dwis_preprocd)
%     f.pick=fnify(dwis_preprocd(i));
    f.pick=dwis_preprocd{i};
    f.no='bvalbvecjsonimg'; %only need to load header
    f.gr=0;
    dwi_clean=d2n2s(folder_in_which_to_write{i},f);
    move_obj_files(dwi_clean,fullfile(folder_in_which_to_write{i},'dke'),working_dwi_name); %I skipped dke, but I still wanna copy the files into a "dke" folder for neatness
    folder_in_which_to_write{i}=fullfile(folder_in_which_to_write{i},'dke'); %after preproc, things should be going in the "final" folder with diffusion metrics and eventually tract outputs
end


end




function [rt1,gmwmi_mask,mask_name,out_5tt]=t1_preproc(t1_folder,dwi_first_b0,out_dir)
%% begin t1 processing
% better to ask for t1 than folder bc t1 metadata isn't required.
%   This goes here instead of beginning bc DWI is loaded now.

% optional: reg mni to t1 before t1 is downsampled?

% note: if the origin of t1 or dwi is whacky, this could cut of some of the
% t1 image. might want to reset origin to center of mass in the future

% load t1
% need to sanitize inputs more
ft1.b0=0; %don't try to assign b0 metadata
ft1.no='bvecbvaljsonfn';
t1=d2n2s(t1_folder,ft1);

% coreg t1 to dwi
% find first b0
% b0_ind_first=paren(find([dwi_dn.bval]<60),1); % get the ind of the first b0

%perform coreg on images in workspace
ft2.apply=1; %only estimate the registration and update the t1's header -- don't reslice the T1
t1=coregister_obj(dwi_first_b0,t1,ft2); % Orientation matrix is changed to align with dwi_dn, but t1 raw data is not changed. (We want full resolution for 5ttgen).


%write
rt1=fullfile(out_dir,'t1_align_with_diff.nii');
[a,b,~]=fileparts(rt1);
d2n2s_write(t1,a,b,[])

%% 5ttgen (and atlas preproc?)

% run 5ttgen
out_5tt=fullfile(out_dir,'5tt.nii'); %I'm not losing anything by using nii and not mif right
system(['5ttgen fsl '...
    rt1 ' '... %in: t1 aligned with diff
    out_5tt ' '... % out: segmentations aligned with diff
    '-nocleanup -force'])

% get a mask output 
temp_dir=clean_dir(dir('./5ttgen-tmp*'));
%grab the most recent folder that matches this pattern
temp_dir=temp_dir(output(@() max(cat(1,temp_dir.datenum)),2)); 

mask_name=[fnify2(temp_dir) filesep 'T1_BET.nii.gz']; %this betted image is output by mrtrix
movefile(mask_name,[out_dir filesep 'T1_BET.nii.gz']) %move the image out of the temp directory
rmdir(fnify2(temp_dir),'s');

% get gmwmi
gmwmi_mask=[out_dir filesep 'gmwmi_mask.nii'];
system(['5tt2gmwmi '...
    out_5tt ' '...
    gmwmi_mask])

% we're done with T1s for now -- we wanted to set up for tractography later, so we
% need to prepare DWIs

end



function []=track_each(folder_in_which_to_write,working_dwi_name,mask_name,gmwmi_mask)
%% get tensors, metrics, SH, and track.
ft_params_template=[fileparts(which(mfilename)) filesep 'ft_parameters.txt']; %just look in the same folder as this function for the template

for i=1:numel(folder_in_which_to_write)

    % define in and out names
    in=[folder_in_which_to_write{i} filesep working_dwi_name];
    in_nii=[in '.nii'];
    in_bval=[in '.bval'];
    in_bvec=[in '.bvec'];
    
    out_DT=[fileparts(in) filesep 'DT.nii'];
    out_KT=[fileparts(in) filesep 'KT.nii'];
    
    %run dwi2tensor
    system(['dwi2tensor ' in_nii ' ' out_DT ' -dkt ' out_KT ' -fslgrad ' in_bvec ' ' in_bval ' -mask ' mask_name])
    
    % tensor2metric or something else for at least fa. probably use
    % designer in the future
    out_FA=[fileparts(in) filesep 'fa.nii'];
    system(['tensor2metric ' out_DT ' -fa ' out_FA  ' -mask ' mask_name])
    
    %kODF for SH_coeff
    study_dirr=fileparts(in);
    out_SH=[folder_in_which_to_write{i} filesep 'dke' filesep 'SH_coeff.nii']; 
    out_tracks=[folder_in_which_to_write{i} filesep 'dke' filesep 'tracks.tck'];
    
    kODF_nii_preprocess3(ft_params_template,out_DT,out_KT,out_FA,[study_dirr filesep],mask_name) % add seed mask but I think this applies only to Euler tracking?
    % the only way in which kODF_nii_preprocess3 is different from the
    % orignal is in that it takes in a studydir which overwrites the ft
    % parameters studydir value. (same with seed mask) Just so i didn't have to write a bunch of
    % ft_parameterses
    

    command=['tckgen -algorithm SD_STREAM -seed_image ' gmwmi_mask ' -mask ' mask_name ' ' out_SH ' ' out_tracks ' -cutoff 0.1 -seeds 10000 -select 10000 -angle 60 -force'];
    system(command)
    
    % just tested this all lightly -- masking with t1 is fine. gmwmi_mask
    % seeding is fine. wm99thresh seeding is fine. they seem similar. but -cutoff 0.05
    % make the tracts look gross.
end


%% generate screenshots of realized tracts
system(['ii=1; for folder in ' out_dir '/sim*/dke; do '...
    'mrview -load $folder/fa.nii '...
    '-tractography.load $folder/tracks.tck '...
    '-capture.folder ' out_dir ' '...
    '-capture.prefix tract_pic$ii '...
    '-capture.grab -exit; '...
    'ii=$(expr $ii + 1); done'])

%% combine .tck files

out_tract_name='all_tracks.tck';
system(['tckedit '...
    out_dir '/sim*/dke/tracks.tck '...
    out_dir filesep out_tract_name]);

%% use tck2connectome to get connectomes
%warp MNI to T1_in_diffusion and atlas along with

% might should make some of these repeated file handling commands functions

%put mni temp in folder
mni_brain_template='/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz'; %gonna have to find this robustly, perhaps using rorden tools
copyfile(mni_brain_template,[out_dir filesep 'mnit1.nii.gz'])
mni_brain_template=[out_dir filesep 'mnit1.nii.gz'];

%gunzip mni temp
try gunzip(mni_brain_template);catch;end
mni_brain_template=strrep(mni_brain_template,'.gz','');

% %gunzip mask % why was I doing this...
% try gunzip(mask_name); catch;end
% mask_name=strrep(mask_name,'.gz','');

%put atlas in folder, gunzip if necessary
[~,b,c]=fileparts(atlas4connectome); %this atlas should be in mni space...
copyfile(atlas4connectome,[out_dir filesep b c]);
atlas4connectome=[out_dir filesep b c];

%gunzip atlas 
try gunzip(atlas4connectome); catch;end
atlas4connectome=strrep(atlas4connectome,'.gz','');

% oldNormSub({mni_brain_template,atlas4connectome},mask_name,8,10,0); %I
% had this line previously... don't know why?
oldNormSub({mni_brain_template,atlas4connectome},rt1,8,10,0); %using nearest neighbor bc labels

%delete moved MNI image
[a,b,c]=fileparts(mni_brain_template);
delete([a filesep 'w' b c])

%change working name of atlas you moved
[a,b,c]=fileparts(atlas4connectome);
atlas4connectome=[a filesep 'w' b c];


out_connectome=[out_dir filesep 'connectome_myc1.csv'];
if exist('atlas4connectome','var') && ~isempty(atlas4connectome) && exist(atlas4connectome,'file')
    system(['tck2connectome '...
        out_dir filesep out_tract_name ' '...
        atlas4connectome ' '...
        out_connectome ' -force'])
end
end


%% LITTLE FUNCTIONS USED

function cdir=clean_dir(dir)
bdir=arrayfun(@(x) x.name(1)=='.',dir);
cdir=dir(~bdir);
end


function oldNormSub(src, tar, smoref, reg, interp)
% usage: oldNormSub(blah1,blah2,8,10,1)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end
if ~exist('smoref','var'), smoref = 0; end
if ~exist('reg','var'), reg = 1; end
if ~exist('interp','var'), interp = 1; end
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src(:);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smoref; % <-- !!!
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = reg;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];%[-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [nan nan nan];%[1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = interp;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
end