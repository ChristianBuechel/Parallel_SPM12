function firstlevel_parallel

hostname = char(getHostName(java.net.InetAddress.getLocalHost));
switch hostname
    case 'mahapralaya'
        base_dir          = 'd:\offhum_cb\';
        n_proc            = 8;
    case 'REM'
        base_dir          = 'c:\Users\buechel\Data\offhum_cb\';
        n_proc            = 2;
    otherwise
        error('Only hosts REM and mahapralaya accepted');
end

%user specified variables
%all_subs    = [9:12 14 15 17:19 21 22 24 26 27 29:35 37:39];
%all_subs  = [9];
%all_subs    = [10:12 14 15 17:19 21 22 24 26 27 29:35 37:39];
%all_subs     = [4:8 16];
all_subs    = [4:12 14 15 16 17:19 21 22 24 26 27 29:35 37:39];

TR                = 1.58;
shift             = 0; %no onset shift

anadirname        = ['first_34_nomov_FAST'];
skern             = 6;

skernel           = repmat(skern,1,3);

struc_templ       = '^sPRISMA.*\.nii';
rfunc_file        = '^ufMRI.nii';


epi_folders       = {'Run1','Run2','Run3'};
conditions        = {'Constant','Offset','Constant_Rate','Offset_Rate'};

n_sess            = size(epi_folders,2);
n_cond            = size(conditions,2);
dummies           = 0;

do_model  = 0;
do_cons   = 0;
do_warp   = 0;
do_smooth = 1;


spm_path          = fileparts(which('spm')); %get spm path
mat_name          = which(mfilename);
[~,mat_name,~]    = fileparts(mat_name);


%prepare for multiprocessing
if size(all_subs) < n_proc
    n_proc = size(all_subs,2);
end 
subs              = splitvect(all_subs, n_proc);

fir_order  = 34;


for np = 1:size(subs,2)
    matlabbatch = [];
    mbi   = 0;
    
    for g = 1:size(subs{np},2)
        %-------------------------------
        %House keeping stuff
        
        name   = sprintf('Sub%02.2d',subs{np}(g));
        name_s = sprintf('sub%02.2d',subs{np}(g));
        st_dir       = [base_dir name filesep 'T1' filesep];
        struc_file   = spm_select('FPList', st_dir, struc_templ);
        for l=1:n_sess
            epi_files{l} = spm_select('ExtFPList', [base_dir filesep name filesep epi_folders{l}], rfunc_file,inf);
        end
        u_rc1_file   = ins_letter(struc_file,'u_rc1');
        
        a_dir    = [base_dir name filesep anadirname];
        
        template = [];
        template.spm.stats.fmri_spec.timing.units   = 'scans';
        template.spm.stats.fmri_spec.timing.RT      = TR;
        template.spm.stats.fmri_spec.timing.fmri_t  = 16;
        template.spm.stats.fmri_spec.timing.fmri_t0 = 1;
        
        template.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        template.spm.stats.fmri_spec.bases.fir.length = TR*fir_order;
        template.spm.stats.fmri_spec.bases.fir.order  = fir_order;
        template.spm.stats.fmri_spec.volt             = 1;
        template.spm.stats.fmri_spec.mthresh          = -Inf;
        template.spm.stats.fmri_spec.global           = 'None';
        %template.spm.stats.fmri_spec.cvi              = 'AR(1)';
        %template.spm.stats.fmri_spec.cvi              = 'None';
        template.spm.stats.fmri_spec.cvi              = 'FAST';
        template.spm.stats.fmri_spec.mask             = cellstr([st_dir 's3skull_strip.nii']);
        
        %template.spm.stats.fmri_spec.cvi              = 'None';
        for sess = 1:n_sess
            s_dir    = [base_dir name filesep epi_folders{sess}];
            fm       = spm_select('FPList', s_dir, '^rp_fMR.*\.txt');
            movement = normit(load(fm));
            %movement = [];
            %all_nuis{sess} = [movement(dummies+1:end,:)];
            all_nuis{sess} = [];
            n_nuis         = size(all_nuis{sess},2);
            z{sess}        = zeros(1,n_nuis); %handy for contrast def
            
            template.spm.stats.fmri_spec.sess(sess).scans = cellstr(epi_files{sess});
            %template.spm.stats.fmri_spec.sess(sess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            template.spm.stats.fmri_spec.sess(sess).multi = {''};
            onset_file = spm_select('FPList',[base_dir 'logfiles' filesep name_s filesep],['^' name_s '_run' num2str(sess) '.*\.mat']);
            RES = extract_offset(onset_file, conditions);
            
            %now do conditions
            for conds = 1:size(RES,2)
                template.spm.stats.fmri_spec.sess(sess).cond(conds).name     = RES{conds}.name;
                template.spm.stats.fmri_spec.sess(sess).cond(conds).onset    = RES{conds}.onset+shift;
                template.spm.stats.fmri_spec.sess(sess).cond(conds).duration = 0;
            end
            
            
            
            template.spm.stats.fmri_spec.sess(sess).multi_reg = {''};
            template.spm.stats.fmri_spec.sess(sess).hpf = 128;
            for nuis = 1:n_nuis
                template.spm.stats.fmri_spec.sess(sess).regress(nuis) = struct('name', cellstr(num2str(nuis)), 'val', all_nuis{sess}(:,nuis));
            end
        end
        
        
        
        if do_model
            mbi = mbi + 1;
            matlabbatch{mbi} = template;
            mkdir(a_dir);
            copyfile(which(mfilename),a_dir);
            matlabbatch{mbi}.spm.stats.fmri_spec.dir = {[a_dir]};
            
            mbi = mbi + 1;
            matlabbatch{mbi}.spm.stats.fmri_est.spmmat           = {[a_dir filesep 'SPM.mat']};
            matlabbatch{mbi}.spm.stats.fmri_est.method.Classical = 1;
        end
                
        
        %%template for contrasts
        template = [];
        template.spm.stats.con.spmmat = {[a_dir filesep 'SPM.mat']};
        template.spm.stats.con.delete = 1;
        fco = 0;
        fco = fco + 1; %counter for f-contrasts
        template.spm.stats.con.consess{fco}.fcon.name   = 'eff_of_int';
        template.spm.stats.con.consess{fco}.fcon.convec = {[repmat([repmat([eye(fir_order)],1,n_cond) zeros(fir_order,n_nuis)],1,n_sess) zeros(fir_order,n_sess)]};
        co_i = 0;
        for co = 1:n_cond
            for i_fir = 1:fir_order
                tpl        = zeros(1,fir_order);
                tpl(i_fir) = 1;
                tpl        = [zeros(1,(co-1)*fir_order) tpl zeros(1,(n_cond-co)*fir_order)];
                convec = [];
                for i_sess = 1:n_sess
                    convec = [convec tpl z{i_sess}];
                end
                co_i = co_i + 1;
                template.spm.stats.con.consess{co_i+fco}.tcon.name    = [conditions{co} '_' num2str(i_fir)];
                template.spm.stats.con.consess{co_i+fco}.tcon.convec  = [convec zeros(1,size(epi_folders,2))];
                template.spm.stats.con.consess{co_i+fco}.tcon.sessrep = 'none';
            end
            
        end
        
        
        
        
        if do_cons
            mbi = mbi + 1;
            matlabbatch{mbi} = template; %now add copnstrasts
        end
        
        
        %prepare_warp
        template = [];
        con_files   = '';
        for co = 1:co_i
            con_files(co,:) = [a_dir filesep sprintf('con_%04.4d.nii',co+fco)];
        end
        
        
        wcon_files             = ins_letter(con_files,'w');
        wcon_dartel_files      = ins_letter(con_files,'w_dartel');
        
        wcon_files             = chng_path(wcon_files, st_dir);    %wcon files still in t1 dir
                
        wcon_dartel_files      = chng_path(wcon_dartel_files, st_dir); %wcon files still in t1 dir
        wcon_dartel_files2     = chng_path(wcon_dartel_files, a_dir);  %wcon files still in ana dir
                
        
        template.spm.tools.dartel.crt_warped.flowfields = cellstr(repmat(u_rc1_file,size(con_files,1),1));
        template.spm.tools.dartel.crt_warped.images = {cellstr(strvcat(con_files))};
        template.spm.tools.dartel.crt_warped.jactransf = 0;
        template.spm.tools.dartel.crt_warped.K = 6;
        template.spm.tools.dartel.crt_warped.interp = 1;
               
        
        if do_warp
            %T1 based dartel
            mbi = mbi + 1;
            matlabbatch{mbi} = template; %now add T1 dartel warp
            
            mbi = mbi + 1;
            matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(wcon_files);
            matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = cellstr(a_dir);
            matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'w';
            matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl    = 'w_dartel';
            matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique         = false;
                                    
        end
        
        if do_smooth
            mbi = mbi + 1;
            matlabbatch{mbi}.spm.spatial.smooth.data = cellstr(wcon_dartel_files2);
            matlabbatch{mbi}.spm.spatial.smooth.fwhm = skernel;
            matlabbatch{mbi}.spm.spatial.smooth.prefix = ['s' num2str(skern)];
        end
    end   
    save([num2str(np) '_' mat_name],'matlabbatch');
    lo_cmd = ['clear matlabbatch;load(''' num2str(np) '_' mat_name ''');'];
    ex_cmd = ['addpath(''' spm_path ''');spm(''defaults'',''FMRI'');spm_jobman(''initcfg'');spm_jobman(''run'',matlabbatch);exit'];
    system(['start matlab.exe -nodesktop -nosplash  -logfile ' num2str(np) '_' mat_name '.log -r "' lo_cmd ex_cmd ';exit"']);
end

function chuckCell = splitvect(v, n)
% Splits a vector into number of n chunks of  the same size (if possible).
% In not possible the chunks are almost of equal size.
%
% based on http://code.activestate.com/recipes/425044/

chuckCell = {};

vectLength = numel(v);


splitsize = 1/n*vectLength;

for i = 1:n
    %newVector(end + 1) =
    idxs = [floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1]+1;
    chuckCell{end + 1} = v(idxs);
end

function out = ins_letter(pscan,letter)
for a=1:size(pscan,1)
    [p , f, e] = fileparts(pscan(a,:));
    out(a,:) = [p filesep letter f e];
end

function out = chng_path(pscan,pa)
for a=1:size(pscan,1)
    [p , f, e] = fileparts(pscan(a,:));
    out(a,:) = [pa filesep f e];
end


function RES = extract_offset(onset_file,conditions)
r1 = load(onset_file);
trials      = r1.p.presentation.trialList;
rate_trials = cell2mat (r1.p.log.onratings.conTrial);
trials(rate_trials) = trials(rate_trials)+2; %make rate trials index 3 and 4
onsets              = r1.p.log.PainOnsetScan
for i=1:size(conditions,2)
    RES{i}.name  = conditions{i};
    RES{i}.onset = onsets(find(trials == i));
end









