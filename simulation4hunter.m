function prob_deterministic_tracking(denoised_dwi_folder,noise_folder,out_dir)

%% the setup

paren=@(x,varargin) x(varargin{:});
fnify=@(x) [x.folder filesep x.name];
curly=@(x,varargin) x{varargin{:}};

% NOTES FOR RUNNING:

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


%read
ff.b0=0;
% dwi_folder='/Users/andrew/re/test/probabilistic_test/denoised'; NOW AN
% INPUT
dwi_dn=d2n2s(dwi_folder,ff);
% noise_folder='/Users/andrew/re/test/probabilistic_test/noise'; NOW AN
% INPUT
noise=d2n2s(noise_folder,ff);

%sizes
nvol=numel(dwi_dn);
dwi_size=size(dwi_dn(1).img);

%constraint and reshape
noise.img(isnan(noise.img))=0;
noise_img=reshape(noise.img,[],1);

%reshape
dwi_dn_img=cat(4,dwi_dn.img);
dwi_dn_img=reshape(dwi_dn_img,[],nvol);


% nSim_total=20;
nsim=10;

% if you wanted to use emilie's tracts
% addpath('/Users/andrew/bin/nii_preprocess2/nii_preprocess/DKI_tractography')

%make a big ole matrix of random numbers
randomm=rand(387200,138,10,2);

%give a working dwi name just to make things a little more programmatic
working_dwi_name='dwi';

%% create and write realizations
tic

noisefunc=@(im_noise,rand_img1,rand_img2) sqrt(2.*im_noise.^2).*(erfinv(2.*rand_img1-1) + 1i.*erfinv(2.*rand_img2-1));
sim_noise_img=noisefunc(noise_img,randomm(:,:,:,1),randomm(:,:,:,2));

signal=abs(dwi_dn_img+sim_noise_img);%signal is now all blowed out in the three d: voxels by volumes by nsims
signal_img = reshape(signal,[dwi_size,nvol,nsim]);
signal_imgc=squeeze(num2cell(signal_img,[1 2 3])); %wonder if there's a faster way...


for u = 1:nsim %this could probably be parallelized -- but passing signal_imgc into the loop would be annoying overhead. Should just make a variable for each segment I think.
    
    fprintf('\n simulation andrew ver #%03d \n',u)

    [dwi_dn.img]=signal_imgc{:,u};
    
    f.gr=0; %this was necessary to put into designer, but you could use f=[] if you wanted
    folder_in_which_to_write{u}=fullfile(out_dir,['sim' num2str(u)]);
    d2n2s_write(dwi_dn,folder_in_which_to_write{u},working_dwi_name,f)    
    
end

reall_time=toc;
fprintf('Time to add noise and write realizations: ',num2str(reall_time))

%% I forgot to do this before, but ya gotta denoise it again
% I can just put it in the folder it was in before and use pick later
appendd=@(x,str) fullfile(choose_output(@() fileparts(x),1),[choose_output(@() fileparts(x),2) str choose_output(@() fileparts(x),3)]); %changed slightly so that you don't need path
preppendd=@(x,str) fullfile(choose_output(@() fileparts(x),1),[str choose_output(@() fileparts(x),2) choose_output(@() fileparts(x),3)]); %changed slightly so that you don't need path

%this is gross rn clean it up later
in=[working_dwi_name '.nii'];
out=appendd(in,'_dn');
cellfun(@(x) system(['dwidenoise ' fullfile(x,in) ' ' fullfile(x,out) ]),folder_in_which_to_write,'un',0)
[~,working_dwi_name,~]=fileparts(out);

%% just put it through DKE
%% not super happy with the file handling going on here, it's just not pretty

% dwis_preprocd=clean_dir(dir('/Users/andrew/re/test/probabilistic_test/sims4dke/*/out/dwi_designer.nii'));
% dwis_preprocd=clean_dir(dir(['/Users/andrew/re/test/probabilistic_test/NEWsim/1sim*andrew/' working_dwi_name '.nii'])); %CHANGE; this should just grab the dwis created by the d2n2s_write command directly above; it matches any instance of folder_in_which_to_write
dwis_preprocd=cellfun(@(x) fullfile(x,out),folder_in_which_to_write,'un',0)'; %this should just grab the dwis created by the d2n2s_write command directly above;
%load into matlab
for i=1:numel(dwis_preprocd)
%     f.pick=fnify(dwis_preprocd(i));
    f.pick=dwis_preprocd{i};
    dwi_clean{i}=d2n2s(folder_in_which_to_write{i},f);
end


in=out;
out=preppendd(in,'dke_');
working_dwi_name=preppendd(working_dwi_name,'dke_'); %just keeping the name separate from the nii.
%make everything pretty the way DKE likes it
for i=1:length(dwi_clean)
    %avg b0s
    b0inds=find([dwi_clean{i}.bval]==0); 
    meanb0img=mean(cat(4,dwi_clean{i}(b0inds).img),4);
    %delete other b0s and replace first with avg
    dwi_clean{i}(b0inds(2:end))=[]; 
    dwi_clean{i}(b0inds(1)).img=meanb0img;
    %remove the b0s bvectors for the grad file that will be made
%     dwi_clean{i}(b0inds(1)).bvec=[]; WE'RE NOT DOING THIS ANY MORE
%     BECAUSE WE'RE NOT USING DKE ANYMORE
    %remove nan and inf
    for j=1:numel(dwi_clean{i})
        dwi_clean{i}(j).img(isnan(dwi_clean{i}(j).img)|isinf(dwi_clean{i}(j).img))=0;
    end
    %write files
    d2n2s_write(dwi_clean{i},[folder_in_which_to_write{i} filesep 'dke'],working_dwi_name,[]) %cut off file extension on out. ; out was previously 
   
end

% in=out;

ft_params_template=['/Users/andrew/re/test/probabilistic_test/las_test/sim1_3/ft_parameters3.txt']; %CHANGE
numb1=find(abs([dwi_clean{i}.bval]-1000)<100);
numb2=find(abs([dwi_clean{i}.bval]-2000)<100);
ndirr=[numel(numb1),numel(numb2)]; 

%should I have a different mask for each realization?
mask_name=[fileparts(folder_in_which_to_write{i}) filesep 'b0.nii']; %
system(['bet ' folder_in_which_to_write{i} filesep 'dke' filesep working_dwi_name ' ' mask_name ' -n -m -f .5'])
mask_name=appendd(mask_name,'_mask'); % bet does this for some reason
system(['fslmaths ' mask_name ' -ero ' mask_name])
mask_name=[mask_name '.gz'];

for i=1:numel(dwis_preprocd)

    % actually just use mrtrix dwi2tensor
    in_nii=[folder_in_which_to_write{i} filesep 'dke' filesep working_dwi_name '.nii'];
    in_bval=[folder_in_which_to_write{i} filesep 'dke' filesep working_dwi_name '.bval'];
    in_bvec=[folder_in_which_to_write{i} filesep 'dke' filesep working_dwi_name '.bvec'];
    
    out_DT=[folder_in_which_to_write{i} filesep 'dke' filesep 'DT' '.nii'];
    out_KT=[folder_in_which_to_write{i} filesep 'dke' filesep 'KT' '.nii'];
    
    system(['dwi2tensor ' in_nii ' ' out_DT ' -dkt ' out_KT ' -fslgrad ' in_bvec ' ' in_bval])
    
    % tensor2metric or something else for at least fa
    out_FA=[folder_in_which_to_write{i} filesep 'dke' filesep 'fa' '.nii'];
    system(['tensor2metric ' out_DT ' -fa ' out_FA])
    
    %kODF for SH_coeff
    study_dirr=[folder_in_which_to_write{i} filesep 'dke'];
    kODF_nii_preprocess3(ft_params_template,out_DT,out_KT,out_FA,[study_dirr filesep],mask_name) % add seed mask but I think this applies only to Euler tracking?
    % the only way in which kODF_nii_preprocess3 is different from the
    % orignal is in that it takes in a studydir which overwrites the ft
    % parameters studydir value. (same with seed mask) Just so i didn't have to write a bunch of
    % ft_parameterses
end


%% mask fibs and do dsi tractography

%grab the fibs
fibs_fn=cellfun(@(x) fullfile(x,'dke','dki.fib') ,folder_in_which_to_write,'un',0);
% use mask_fib on all
cellfun(@(x) mask_fib(mask_name,x,[fileparts(x) filesep 'masked.fib']),fibs_fn)

dsi_location='/Applications/dsi_studio.app/Contents/MacOS/dsi_studio'; %CHANGE
for i=1:length(fibs_fn)
%     system([dsi_location ' --action=trk --source=' fibs_dir(i).folder filesep 'masked.fib --connectivity=aal --connectivity_value=count,ncount,ncount2 --seed_count=1000000 --fa_threshold=0.1 --turning_angle=35 --step_size=1 --smoothing=0 --connectivity_type=pass,end --output=' fibs_dir(i).folder filesep 'track_results.txt --export=stat'])
    system([dsi_location ' --action=trk --source=' fileparts(fibs_fn{i}) filesep 'masked.fib --connectivity=FreeSurferSeg --connectivity_value=count --seed_count=1000000 --fa_threshold=0.1 --turning_angle=35 --step_size=1 --smoothing=0 --min_length=30 --seed_plan=1 --initial_dir=2 --connectivity_type=pass,end --output=' fileparts(fibs_fn{i}) filesep 'track_results.txt --export=stat'])

end

for i=1:10
%     system([dsi_location ' --action=trk --source=' fibs_dir(i).folder filesep 'masked.fib --connectivity=aal --connectivity_value=count,ncount,ncount2 --seed_count=1000000 --fa_threshold=0.1 --turning_angle=35 --step_size=1 --smoothing=0 --connectivity_type=pass,end --output=' fibs_dir(i).folder filesep 'track_results.txt --export=stat'])
    fprintf([dsi_location ' --action=trk --source=' fileparts(fibs_fn{1}) filesep 'masked.fib --connectivity=FreeSurferSeg --connectivity_value=count,ncount,ncount2 --seed_count=1000000 --fa_threshold=0.1 --turning_angle=35 --step_size=1 --smoothing=0 --min_length=30 --seed_plan=1 --initial_dir=2 --connectivity_type=pass,end --output=' fileparts(fibs_fn{1}) filesep 'simSIM' filesep 'track_results' num2str(i) '.txt --export=stat\n'])
    
end


%% random stuff related to connectomes; just messing around

% cgs_fn=clean_dir(dir('/Users/andrew/re/test/probabilistic_test/sims4dke/*/out/dke/track_results.txt.aal.count.pass.connectogram.txt'));
% connectome=load_connectograms('/Users/andrew/re/test/probabilistic_test/sims4dke/*/out/dke/track_results.txt.aal.count.pass.connectogram.txt');
connectome=load_connectograms([fileparts(folder_in_which_to_write{1}) filesep '*/dke/track_results.txt.*.count.pass.connectogram.txt']);
connectome=cat(3,connectome{:});

% imagesc(cell2mat(connectome{1}(3:end,3:end)))
% con_mat=cellfun(@(x) cell2mat(x(3:end,3:end)),connectome,'un',0);
ct_mat=cell2mat(connectome(3:end,3:end,:));
figure
for i=1:size(ct_mat,3)
%     figure
    imagesc(log(ct_mat(:,:,i)))
    pause(3)
end
%that actually shows some variability

ct_mean=mean(ct_mat,3);
ct_sd=std(ct_mat,0,3);

figure
imagesc(ct_mean)
imagesc(ct_sd)

median(ct_sd(ct_mean~=0))
median(ct_mean(ct_mean~=0))

% 
% 
% ct_mat()
% 
% isnormall=cellfun(@(x) adtest(squeeze(x)),num2cell(ct_mat,3));
% ct_vecs=num2cell(ct_mat,3);
% for i=1:numel(ct_vecs)
%     testss=adtest(ct_vecs{i});
% end




%% other functions are below. some suck, don't judge

function mask_fib(mask_fn,fib_fn,out_fn)
% ugly script, but it works. and the methods used are
% similar to those suggested by frank yeh of DSI studio on the website ...
% or just copied from it lol
%


gunziped=@(x) x(1:end-3);
paren=@(x,varargin) x(varargin{:});

if contains(mask_fn,'.gz')
    gunzip(mask_fn)
    mask_fn=gunziped(mask_fn);
end

%load
mask_hdr=spm_vol(mask_fn);
mask_img=spm_read_vols(mask_hdr(1));
if size(mask_img,4)~=1; warning('your input mask is 4d -- function is intended to be used with 3d data');end
mask_img(mask_img~=0)=1; % binarize. this is idempotent


%Image must be LPS -- that's all DSI studio takes, apparently. So we
%probably need to flip the mask images 

%flip a/p if needed
if mask_hdr(1).mat(6)>0
    mask_img=mask_img(:,size(mask_img,2):-1:1,:);
end
%flip l/r if needed
if mask_hdr(1).mat(1)>0
    mask_img=mask_img(size(mask_img,1):-1:1,:,:);
end

%flip i/s if needed
if mask_hdr(1).mat(11)<0
    mask_img=mask_img(:,:,size(mask_img,1):-1:1);
end

fib = load(fib_fn,'-mat');

%find number of each var
max_fib = 0;
for i = 1:10
    if isfield(fib,strcat('fa',int2str(i-1))) 
        max_fib = i;
    else
        break;
    end
end

max_odf = 0;
for i = 1:10
    if isfield(fib,strcat('odf',int2str(i-1))) 
        max_odf = i;
    else
        break;
    end
end




nn=size(fib.odf0,2); %size of chunks

fa_nz_vector=fib.fa0~=0;
fa_nz_linear_indices=find(fib.fa0~=0); %inds of original
mask_img_vector=reshape(mask_img,1,[]);

deletion_vector=fa_nz_vector & ~mask_img_vector;
deletion_linear_indices=find(deletion_vector); %"missing" inds
disp(['numel deleted = ' num2str(sum(deletion_vector(:)))])

%do the magic
logical_for_odf_deletion=ismember(fa_nz_linear_indices,deletion_linear_indices);


odf=[];
for i = 1:max_odf
    %build odf var
    eval(strcat('odf =  cat(2,odf,fib.odf',int2str(i-1),');'))    
    strcat('odf =  cat(2,odf,fib.odf',int2str(i-1),');')

end

if size(odf,2)~=sum(fib.fa0(:)~=0)
    warning('odf and fa number of elements is not the same')
end

odf(:,logical_for_odf_deletion)=[]; 
for i = 1:max_fib
    eval(strcat('fa',int2str(i-1),'=fib.fa',int2str(i-1),';'));
    strcat('fa',int2str(i-1),'=fib.fa',int2str(i-1),';')
    eval(strcat('dir',int2str(i-1),'=fib.dir',int2str(i-1),';'));
    strcat('dir',int2str(i-1),'=fib.dir',int2str(i-1),';')
    eval(strcat('fa',int2str(i-1),'(1,deletion_vector)=0;',';'))
    strcat('fa',int2str(i-1),'(1,deletion_vector)=0;')
    eval(strcat('dir',int2str(i-1),'(:,deletion_vector)=0;',';'))
    strcat('dir',int2str(i-1),'(:,deletion_vector)=0;')
%     eval(strcat('clear odf',int2str(i-1))) % unnecessary -- so far, we
%     haven't defined an odf_i outside of fib.odf_i
%     strcat('clear odf',int2str(i-1))
end

%rewrite odf
%how many chunks do we need?
num_odf_vars_needed=ceil(numel(odf)/(nn*size(fib.odf0,1))); 


for i=1:num_odf_vars_needed
    starting_index=(i-1).*nn+1;
    ending_index=(i).*nn;
    if i~=num_odf_vars_needed
        eval(strcat('odf',int2str(i-1),'=odf(:,starting_index:ending_index);'));
        strcat('odf',int2str(i-1),'=odf(:,starting_index:ending_index);')
    elseif i==num_odf_vars_needed
        eval(strcat('odf',int2str(i-1),'=odf(:,starting_index:end);')); %could just find the remainder instead of doing 'end', but whatever
        strcat('odf',int2str(i-1),'=odf(:,starting_index:end);')
    end
end
%keep only the variables you want to save
odf_faces=fib.odf_faces;
odf_vertices=fib.odf_vertices;
dimension=fib.dimension;
voxel_size=fib.voxel_size;

clear odf
clear mask_img_vector
clear fa_nz_vector
clear deletion_vector
clear mask_img
clear max_fib
clear i
clear starting_index % surely there's a better way than this
clear ending_index   % you could just do multiple save(...,'-append') statements wrapped in eval(...) to refer to the variable number of odf,fa, etc variables
clear fib
clear deletion_vector
clear deletion_linear_indices
clear fa_nz_linear_indices
clear logical_for_odf_deletion
clear num_odf_vars_needed
clear nn
clear gunziped
clear mask_fn
clear paren
clear ans
clear fib_fn
clear max_odf
clear mask_hdr


save(out_fn,'-v4') %this is going to save out_fn tho.... might be okay.
end

function dke_options(varargin)

% dke Diffusional Kurtosis Estimator
%   inputs are dicom or nifti diffusion-weighted images 
%   outputs are diffusion and diffusional kurtosis tensors and tensor-derived maps

% Author: Ali Tabesh
% Version: 2.6.0
% Last modified: 11/5/2014 by EM

% -------------------------------------------
% Constants
% -------------------------------------------

% if set to 1, dke will output tensors, e'vecs, and e'vals, 
% as well as d2/k2, d3/k3, color fa, and white matter model-derived maps
internal_flag = 1;

% -------------------------------------------
% Check inputs
% -------------------------------------------

[dkeVersion, dkeDate] = GetDKEVersion;
fprintf('%s, %s\n', dkeVersion, dkeDate)

% 
% if nargin ~= 1
%     fprintf('\n')
%     fprintf('Usage: dke paramsfile\n')
%     fprintf('paramsfile  DKE processing parameters file\n\n')
%     return
% end
% 
% if ~exist(fn_params, 'file')
%     fprintf('\n')
%     fprintf('Input parameters file %s does not exist!\n\n', fn_params)
%     return
% end

% -------------------------------------------
% Read parameters file
% -------------------------------------------
% 
% fid=fopen(fn_params); %EM
% file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', ''); %EM
% for i = 1:length(file{1}) %EM
%     eval(file{1}{i})%EM
% end
% fclose(fid);
% 
% if exist('subject_list','var') %EM
% else
% subject_list={''};
% end

%-----------------------------
% overwrite parameters that you don't like using input parser (AT,
% 1/29/2019)
%-----------------------------
p=inputParser;

addParameter(p,'studydir','')
addParameter(p,'subject_list',{''})
addParameter(p,'preprocess_options_format','nifti')
addParameter(p,'preprocess_options_navg','1')
addParameter(p,'preprocess_options_extra_b0','1')
addParameter(p,'preprocess_options_coreg_flag','1')
addParameter(p,'preprocess_options_series_description',{'DKI'})
addParameter(p,'preprocess_options_fn_nii','dkis.nii')
addParameter(p,'bval',[0 1000 2000])
addParameter(p,'ndir',[128,128])
addParameter(p,'fn_gradients','')
addParameter(p,'idx_gradients',{1:128;1:128})
addParameter(p,'idx_1st_img',1)
addParameter(p,'Kmin',0)
addParameter(p,'NKmax',3)
addParameter(p,'Kmin_final',0)
addParameter(p,'Kmax_final',3)
addParameter(p,'T',500)
addParameter(p,'find_brain_mask_flag',1)
addParameter(p,'dki_method_no_tensor',0)
addParameter(p,'dki_method_linear_weighting',1)
addParameter(p,'dki_method_linear_constrained',1)
addParameter(p,'dki_method_nonlinear',0)
addParameter(p,'dki_method_linear_violations',0)
addParameter(p,'dki_method_robust_option',0)
addParameter(p,'dki_method_noise_tolerance',0.09)
addParameter(p,'dti_method_dti_flag',0)
addParameter(p,'dti_method_dti_only',0)
addParameter(p,'dti_method_no_tensor',0)
addParameter(p,'dti_method_linear_weighting',1)
addParameter(p,'dti_method_b_value',1000)
addParameter(p,'dti_method_directions',{1:128})
addParameter(p,'dti_method_robust_option',0)
addParameter(p,'dti_method_noise_tolerance',0.09)
addParameter(p,'fwhm_img',[0 0 0])
addParameter(p,'fn_noise','')
addParameter(p,'fwhm_noise',0)
addParameter(p,'median_filter_method',1)
addParameter(p,'map_interpolation_method_flag',0)
addParameter(p,'map_interpolation_method_order',1)
addParameter(p,'map_interpolation_method_resolution',1)

% THIS IS WHERE add....s go!

parse(p,varargin{:});

% toplevelfield='p.Results';
% passed_value_names=fieldnames(p.Results);
% potential_fldnms={'studydir';'subject_list';'preprocess_option.format';'preprocess_options.fn_nii';'bval';'ndir';'fn_gradients';'idx_gradients';'T';'find_brain_mask_flag';'fwhm_img'};

% parsed_fieldnames=cellfun(@(x) [toplevelfield '.' x],potential_fldnms,'UniformOutput',0);

potential_fldnms={'studydir','subject_list','preprocess_options.format','preprocess_options.navg','preprocess_options.extra_b0',' preprocess_options.coreg_flag','  preprocess_options.series_description',' preprocess_options.fn_nii','bval','ndir','fn_gradients','idx_gradients','idx_1st_img','Kmin','NKmax','Kmin_final','Kmax_final','T','find_brain_mask_flag','dki_method.no_tensor','dki_method.linear_weighting','dki_method.linear_constrained','dki_method.nonlinear','dki_method.linear_violations','dki_method.robust_option','dki_method.noise_tolerance','dti_method.dti_flag','dti_method.dti_only','dti_method.no_tensor','dti_method.linear_weighting','dti_method.b_value','dti_method.directions','dti_method.robust_option','dti_method.noise_tolerance','fwhm_img','fn_noise','fwhm_noise','median_filter_method','map_interpolation_method.flag','map_interpolation_method.order',' map_interpolation_method.resolution'};

for i=1:numel(potential_fldnms) %next three lines come from stackexchange
    %parse parsed_fieldnames into references to the actual variables
    if ~contains(potential_fldnms{i},'.')
        eval(sprintf('%s=%s',potential_fldnms{i},['p.Results.' potential_fldnms{i}]))
    elseif contains(potential_fldnms{i},'.')%numel(curly(textscan(potential_fldnms,'%s','Delimiter','.'),1))
        fieldss=textscan(potential_fldnms{i},'%s','Delimiter','.');
        eval(sprintf('%s=%s',potential_fldnms{i},['p.Results.' fieldss{1}{1} '_' fieldss{1}{2}] ))
    end
end

save([studydir filesep 'dke_params.mat'])

% -------------------------------------------
% Output file names
% -------------------------------------------

fn_out_struc.dtype_out = 'float32';                 % output data format

% DKI outputs
fn_out_struc.kmean = 'kmean';                       % Mean kurtosis
fn_out_struc.k2    = 'k2';                          % Kurtosis along direction of medium diffusion
fn_out_struc.k3    = 'k3';                          % Kurtosis along direction of minimum diffusion
fn_out_struc.kax   = 'kax';                         % Axial kurtosis
fn_out_struc.krad  = 'krad';                        % Radial kurtosis

fn_out_struc.dmean        = 'dmean';                % Mean diffusivity
fn_out_struc.d2           = 'd2';                   % Diffusivity along direction of medium diffusion
fn_out_struc.d3           = 'd3';                   % diffusivity along direction of minimum diffusion
fn_out_struc.dax          = 'dax';                  % axial diffusivity
fn_out_struc.drad         = 'drad';                 % radial diffusivity
fn_out_struc.fa           = 'fa';                   % fractional anisotropy (FA)
fn_out_struc.fa_color_map = 'fa_color_map';         % FA color map

fn_out_struc.nexcl           = 'noutlier';          % number of outliers with the robust method (with method.robust_option ~= 0)
fn_out_struc.fit_err         = 'fit_err';           % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_avg     = 'fit_err_avg';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_map = 'abs_fit_err_map';   % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_map     = 'fit_err_map';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_avg = 'abs_fit_err_avg';   % fraction of unexplained signal variance by DKI model

fn_out_struc.d_viol    = 'd_viol';                  % fraction of constraint violations on directional diffusivity
fn_out_struc.kmin_viol = 'kmin_viol';               % fraction of constraint violations on minimum directional kurtosis
fn_out_struc.kmax_viol = 'kmax_viol';               % fraction of constraint violations on maximum directional kurtosis

fn_out_struc.kt  = 'KT.mat';                        % kurtosis tensor
fn_out_struc.dt  = 'DT.mat';                        % DKI diffusion tensor
fn_out_struc.kfa = 'kfa';                           % Kurtosis FA
fn_out_struc.meankt='mkt';                          % Mean Kurtosis Tensor 

fn_out_struc.dt_eigval = 'dt_eigval';               % DKI diffusion tensor eigenvalues
fn_out_struc.dt_eigvec = 'dt_eigvec';               % DKI diffusion tensor eigenvectors

% DTI outputs
fn_out_struc.dmean_dti = 'dmean_dti';               % mean diffusivity from DTI computation
fn_out_struc.d2_dti    = 'd2_dti';                  % diffusivity along direction of medium diffusion from DTI computation
fn_out_struc.d3_dti    = 'd3_dti';                  % diffusivity along direction of minimum diffusion from DTI computation
fn_out_struc.dax_dti   = 'dax_dti';                 % axial diffusivity from DTI computation
fn_out_struc.drad_dti  = 'drad_dti';                % radial diffusivity from DTI computation
fn_out_struc.fa_dti    = 'fa_dti';                  % FA from DTI computation

fn_out_struc.fa_color_map_dti    = 'fa_color_map_dti';      % FA color map
fn_out_struc.nexcl_dti           = 'noutlier_dti';          % number of outliers with the robust method (with method.robust_option ~= 0) from DTI computation
fn_out_struc.fit_err_dti         = 'fit_err_dti';           % fraction of unexplained signal variance by DTI model
fn_out_struc.fit_err_avg_dti     = 'fit_err_avg_dti';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_map_dti = 'abs_fit_err_map_dti';   % fraction of unexplained signal variance by DKI model
fn_out_struc.fit_err_map_dti     = 'fit_err_map_dti';       % fraction of unexplained signal variance by DKI model
fn_out_struc.abs_fit_err_avg_dti = 'abs_fit_err_avg_dti';   % fraction of unexplained signal variance by DKI model

fn_out_struc.dt_dti = 'dt_dti';                     % DTI diffusion tensor
fn_out_struc.dt_dti_eigval = 'dt_dti_eigval';       % DTI diffusion tensor eigenvalues
fn_out_struc.dt_dti_eigvec = 'dt_dti_eigvec';       % DTI diffusion tensor eigenvectors

% WMM outputs

fn_out_struc.awf = 'wmm_awf';                       % WMM axonal water fraction
fn_out_struc.da = 'wmm_da';                         % WMM axonal diffusivity
fn_out_struc.de_axial = 'wmm_de_ax';                % WMM axial extra-axonal diffusivity
fn_out_struc.de_radial = 'wmm_de_rad';              % WMM radial extra-axonal diffusivity
fn_out_struc.tortuosity = 'wmm_tort';               % WMM tortuosity

% -------------------------------------------
% Turn all warnings off
% -------------------------------------------

 warning('off','all');

% -------------------------------------------
% Process all subjects
% -------------------------------------------

for isubject = 1:length(subject_list)

    dir_subj = fullfile(studydir, subject_list{isubject});      % subject root folder

    diary off
    fn_diary = fullfile(dir_subj, 'dke.log');
    if exist(fn_diary, 'file')
        delete(fn_diary);
    end
    fid = fopen(fn_diary, 'w');
    if fid < 0
        error('Cannot open output file %s! Output directory does not exist or is write-protected.', fn_diary);
    end
    diary(fn_diary)

    fprintf('Start date and time: %s\n', datestr(now, 'mmmm dd, yyyy HH:MM:SS'))
    fprintf('%s\n',dkeVersion)% EM
    fnames = fieldnames(fn_out_struc);
    if isempty(subject_list{isubject})
        for ifield = 1:length(fnames)
            eval(['fn_out_struc_subject.' fnames{ifield} ' = fullfile(dir_subj, fn_out_struc.' fnames{ifield} ');']);
        end
    else
        for ifield = 1:length(fnames)
            eval(['fn_out_struc_subject.' fnames{ifield} ' = fullfile(dir_subj, [subject_list{isubject} ''_'' fn_out_struc.' fnames{ifield} ']);']);
        end
    end
    fn_out_struc_subject.dtype_out = fn_out_struc.dtype_out;
    
    if ~isempty(fn_noise)
        fn_subject_noise = fullfile(dir_subj, fn_noise);
    else
        fn_subject_noise = '';
    end
    
    if strcmpi(preprocess_options.format, 'dicom')
        eval(['!dke_preprocess_dicom "' dir_subj '" "' fn_params '"']);
        fn_subject_img_prefix = fullfile(dir_subj, 'intermediate_processing', 'combined', 'rdki');
        
    elseif strcmpi(preprocess_options.format, 'nifti')
        fn_subject_img_prefix = fullfile(dir_subj, preprocess_options.fn_nii);
    elseif strcmpi(preprocess_options.format, 'bruker')
        dke_preprocess_bruker(dir_subj, bval, preprocess_options);
        fn_subject_img_prefix = fullfile(dir_subj, 'intermediate_processing', 'combined', 'rdki');
    else
        error('Invalid input image format! Supported formats are DICOM and NIfTI.')
    end
    
    dke_estimate(fn_subject_img_prefix, idx_1st_img, bval, ndir, Kmin, NKmax, Kmin_final, Kmax_final, ...
        T, find_brain_mask_flag, dki_method, dti_method, fwhm_img, fn_subject_noise, fwhm_noise, ...
        fn_gradients, idx_gradients, fn_out_struc_subject, median_filter_method, map_interpolation_method, internal_flag)

end

diary off
end

function cdir=clean_dir(dir)
bdir=arrayfun(@(x) x.name(1)=='.',dir);
cdir=dir(~bdir);
end

function varargout=load_connectograms(fn_pattern)

% If given, varargout is the results structure array, and the user can make
% sure the index of the filenames matches the index of the matrix cell they
% want to look at 

% Writing this script is just kinda practice with the variable arg tools

nargoutchk(1,2)

%results=dir('V:\taylor_j\leo\dr_cooper\PD_MRIs2\PD_MRIs\VaP*\dke\VaP*_results.txt.atlas_VaP*.count.end.connectogram.txt');
results=dir(fn_pattern);

%sanitize input 'r whatever
problem_index=[];
problem_index2=[];
for j=1:length(results)
    if results(j).isdir==1
        problem_index=[problem_index j];
    end
    if ~contains(results(j).name,'connectogram')
        problem_index2=[problem_index2 j];
    end
end
problem_indices=union(problem_index,problem_index2);
results(problem_indices)=[];

for i=1:length(results)
    con_matrix{i}=delimread([results(i).folder filesep results(i).name],{' ','\t'},'mixed');
    matrix{i}=con_matrix{i}.mixed;
end





varargout{1}=matrix;

if nargout==2
    varargout{2}=results;
end
end


% 
% %% test
% %% does paralellized writing help?
% %% not
% tic
% for u = 1:nsim %this could probably be parallelized -- but passing signal_imgc into the loop would be annoying overhead. Should just make a variable for each segment I think.
%     
%     fprintf('\n simulation andrew ver #%03d \n',u)
% 
%     [dwi_dn.img]=signal_imgc{:,u};
%     
%     f=struct();
%     f.gr=0; %this was necessary to put into designer, but you could use f=[] if you wanted
%     folder_in_which_to_write=['/Users/andrew/re/test/probabilistic_test/NEWsim/1sim' num2str(u) 'andrew4delete']; %CHANGE;
%     d2n2s_write(dwi_dn,folder_in_which_to_write,working_dwi_name,f)    
%     
% end
% fprintf('NOT finished')
% toc
% %% is
% tic
% for u=1:nsim
%     dwi_dn2{u}=dwi_dn;
%     [dwi_dn2{u}.img]=signal_imgc{:,u};
% end
% parfor u = 1:nsim %this could probably be parallelized -- but passing signal_imgc into the loop would be annoying overhead. Should just make a variable for each segment I think.
%     
%     fprintf('\n simulation andrew ver #%03d \n',u)
% 
%     [dwi_dn2{u}.img]=signal_imgc{:,u};
%     
%     f=struct();
%     f.gr=0; %this was necessary to put into designer, but you could use f=[] if you wanted
%     folder_in_which_to_write=['/Users/andrew/re/test/probabilistic_test/NEWsim/1sim' num2str(u) 'andrew4delete2']; %CHANGE;
%     d2n2s_write(dwi_dn2{u},folder_in_which_to_write,working_dwi_name,f)    
%     
% end
% fprintf('IS finished')
% toc