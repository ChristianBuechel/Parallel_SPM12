function preprocessing_fmri

hostname = char(getHostName(java.net.InetAddress.getLocalHost));
switch hostname
    case 'mahapralaya'
        base_dir          = 'd:\offhum_cb\';
        n_proc            = 4;
    case 'REM'
        base_dir          = 'c:\Users\buechel\Data\offhum_cb\';
        n_proc            = 2;
    otherwise
        error('Only hosts REM and mahapralaya accepted');
end

check          = 0;

do_4d          = 0;
do_field       = 0;
do_realign     = 0;
do_real_unwarp = 0;
do_coreg       = 0;
do_seg         = 0;
do_skull       = 0;
do_sm_skull    = 0;
do_norm        = 0;
do_back        = 0;
do_warp        = 0;
do_avg_norm    = 1;

all_subs    = [9:12 14 15 17:19 21 22 24 26 27 29:35 37:39];
%all_subs     = 9:12;
all_subs     = [4:8 16];

all_subs    = [4:8 9:10 12 14 15 16 17:19 21  24 26 27 29:35 37:39]; %exclude 11 22

TR           = 1.58;

sfunc_templ       = '^fPRISMA.*\.nii';
struc_templ       = '^sPRISMA.*\.nii';

%fm_magn_1_templ   = '^sPRISMA.*\-0012.*\-01.*\.nii';
%fm_magn_2_templ   = '^sPRISMA.*\-0012.*\-02.*\.nii';
%fm_phase_templ    = '^sPRISMA.*\-0013.*\-02.*\.nii';
%Volunteers 27 and 31 have -0013 and -0014

func_name         = 'fMRI.nii';

mean_func_name    = 'meanufMRI.nii';

skullstrip_name   = 'skull_strip.nii';

fm_default_file   = 'pm_defaults_Prisma_ISN_15.m';
epi_folders       = {'Run1','Run2','Run3'};
dummies           = 0; %already taken care of


% h=clock;
% while h(4)<18 %wait until 18h to start this
%     h=clock;
% end

if size(all_subs) < n_proc
    n_proc = size(all_subs,2);
end

subs     = splitvect(all_subs, n_proc);
spm_path = fileparts(which('spm')); %get spm path
template_path = [spm_path filesep 'toolbox\cat12\templates_1.50mm' filesep];
tpm_path = [spm_path filesep 'tpm' filesep];


for np = 1:size(subs,2)
    matlabbatch = [];
    gi   = 1;
    
    for g = 1:size(subs{np},2)
        %-------------------------------
        %House keeping stuff
        
        name = sprintf('Sub%02.2d',subs{np}(g));
        st_dir       = [base_dir name filesep 'T1' filesep];
        fm_dir       = [base_dir name filesep 'FM' filesep];
        struc_file   = spm_select('FPList', st_dir, struc_templ);
        
        %get all images in FM
        fm_files     = spm_select('FPList', fm_dir, struc_templ);
        fm_hdr       = spm_data_hdr_read(fm_files);
        min_scale    = cat(2,fm_hdr.pinfo);
        phase_ind    = find(min_scale(2,:)<0); %phase image has negative scale factor at least now ...
        mag_ind      = find(min_scale(2,:)==0); %phase image has negative scale factor at least now ...
        phase_file   = fm_files(phase_ind,:);
        
        %now get magnitude images
        mag_info       = dis_fname(fm_files(mag_ind,:));
        mag_info_count = cat(2,mag_info.count);
        [small_mag,~] = min(mag_info_count);
        [big_mag,~]   = max(mag_info_count);
        
        magn_1_file   = fm_files(small_mag,:);
        magn_2_file   = fm_files(big_mag,:);
        
        for l=1:size(epi_folders,2)
            epi_files{l} = spm_select('FPList', [base_dir filesep name filesep epi_folders{l}], sfunc_templ);
        end
        
        
        if check
            fprintf(['***Volunteer ' name '\n']);
            
            fprintf('Structural image: ');
            if isempty(struc_file)
                fprintf('missing\n')
            else
                fprintf('OK %d file(s)\n',size(struc_file,1))
            end
            
            fprintf('Phase image     : ');
            if isempty(phase_file)
                fprintf('missing\n')
            else
                fprintf('OK %d file(s)\n',size(phase_file,1))
            end
            
            fprintf('Magn image      : ');
            if isempty(magn_1_file)
                fprintf('missing\n')
            else
                fprintf('OK %d file(s)\n',size(magn_1_file,1))
            end
            
            fprintf('EPI images\n');
            for l=1:size(epi_folders,2)
                fprintf('S%d  : ',l);
                if isempty(epi_files{l})
                    fprintf('missing\n');
                else
                    fprintf('OK: %d volumes\n',size(epi_files{1},1));
                end
            end
            fprintf('\n')
            
        end
        
        m_dir        = [base_dir name filesep epi_folders{1} ];
        mean_file    = [m_dir filesep mean_func_name];
        c1_file      = ins_letter(struc_file,'c1');
        c2_file      = ins_letter(struc_file,'c2');
        rc1_file     = ins_letter(struc_file,'rc1');
        rc2_file     = ins_letter(struc_file,'rc2');
        u_rc1_file   = ins_letter(struc_file,'u_rc1');
        strip_file   = [base_dir name filesep 'T1' filesep skullstrip_name];
        
        fprintf(['Doing volunteer ' name '\n']);
        
        %-------------------------------
        %Do 4D NIFTI conversion
        if do_4d
            for l=1:size(epi_folders,2)
                matlabbatch{gi}.spm.util.cat.vols  = cellstr(epi_files{l}(dummies+1:end,:));
                matlabbatch{gi}.spm.util.cat.name  = func_name;
                matlabbatch{gi}.spm.util.cat.dtype = 0;
                matlabbatch{gi}.spm.util.cat.RT    = TR;
                gi = gi + 1;
                %                 matlabbatch{gi}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(epi_files);
                %                 matlabbatch{gi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
                %                 gi = gi + 1;
            end
        end
        %save('4D.mat','matlabbatch');
        
        
        %-------------------------------
        %Do Fieldmap
        if do_field
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {phase_file};
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {magn_1_file};
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile = {fm_default_file};
            for l=1:size(epi_folders,2)
                s_dir            = [base_dir name filesep epi_folders{l}];
                matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.session(l).epi = create_func_files(s_dir,func_name,1)';
            end
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
            matlabbatch{gi}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
            gi = gi + 1;
        end
        
        %-------------------------------
        %Do Realignment
        if do_realign
            
            for l=1:size(epi_files,2)
                s_dir            = [base_dir name filesep epi_folders{l}];
                func_files{l}    = create_func_files(s_dir,func_name,size(epi_files{l},1))';
            end
            
            matlabbatch{gi}.spm.spatial.realign.estwrite.data             = func_files;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
            matlabbatch{gi}.spm.spatial.realign.estwrite.eoptions.weight  = '';
            matlabbatch{gi}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
            matlabbatch{gi}.spm.spatial.realign.estwrite.roptions.interp  = 4;
            matlabbatch{gi}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
            matlabbatch{gi}.spm.spatial.realign.estwrite.roptions.mask    = 1;
            matlabbatch{gi}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
            gi = gi + 1;
        end
        %save('realign.mat','matlabbatch');
        
        %-------------------------------
        %Do Realignment and unwarp
        if do_real_unwarp
            for l=1:size(epi_folders,2)
                s_dir            = [base_dir name filesep epi_folders{l}];
                matlabbatch{gi}.spm.spatial.realignunwarp.data(l).scans = create_func_files(s_dir,func_name,size(epi_files{l},1))';
                matlabbatch{gi}.spm.spatial.realignunwarp.data(l).pmscan = {ins_letter(phase_file,'vdm5_sc',['_session' num2str(l)])};
            end
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.sep = 4;
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.rtm = 0;
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.einterp = 2;
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
            matlabbatch{gi}.spm.spatial.realignunwarp.eoptions.weight = '';
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.jm = 0;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.sot = [];
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.rem = 1;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.noi = 5;
            matlabbatch{gi}.spm.spatial.realignunwarp.uweoptions.expround = 'First';
            matlabbatch{gi}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
            matlabbatch{gi}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
            matlabbatch{gi}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
            matlabbatch{gi}.spm.spatial.realignunwarp.uwroptions.mask = 1;
            matlabbatch{gi}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
            gi = gi + 1;
            
        end
        %-------------------------------
        %Do Coregistration mean rfMRI to T1
        if do_coreg
            matlabbatch{gi}.spm.spatial.coreg.estimate.source = cellstr(struc_file);
            matlabbatch{gi}.spm.spatial.coreg.estimate.ref    = cellstr(mean_file);
            matlabbatch{gi}.spm.spatial.coreg.estimate.other  = {''};
            matlabbatch{gi}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{gi}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{gi}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{gi}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            gi = gi + 1;
        end
        %-------------------------------
        %Do Segmentation
        if do_seg
            matlabbatch{gi}.spm.spatial.preproc.channel.vols     = cellstr(struc_file);
            matlabbatch{gi}.spm.spatial.preproc.channel.biasreg  = 0.001;
            matlabbatch{gi}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{gi}.spm.spatial.preproc.channel.write    = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(1).tpm    = {[tpm_path 'enhanced_TPM.nii,1']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(1).ngaus  = 2; %LOrio ... Draganski et al. NI2016
            matlabbatch{gi}.spm.spatial.preproc.tissue(1).native = [1 1];
            matlabbatch{gi}.spm.spatial.preproc.tissue(1).warped = [0 0];   
            matlabbatch{gi}.spm.spatial.preproc.tissue(2).tpm    = {[tpm_path 'enhanced_TPM.nii,2']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(2).ngaus  = 1;
            matlabbatch{gi}.spm.spatial.preproc.tissue(2).native = [1 1];
            matlabbatch{gi}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(3).tpm    = {[tpm_path 'enhanced_TPM.nii,3']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(3).ngaus  = 2;
            matlabbatch{gi}.spm.spatial.preproc.tissue(3).native = [1 1];
            matlabbatch{gi}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(4).tpm    = {[tpm_path 'enhanced_TPM.nii,4']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(4).ngaus  = 3;
            matlabbatch{gi}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(5).tpm    = {[tpm_path 'enhanced_TPM.nii,5']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(5).ngaus  = 4;
            matlabbatch{gi}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(6).tpm    = {[tpm_path 'enhanced_TPM.nii,6']};
            matlabbatch{gi}.spm.spatial.preproc.tissue(6).ngaus  = 2;
            matlabbatch{gi}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{gi}.spm.spatial.preproc.warp.mrf         = 1;
            matlabbatch{gi}.spm.spatial.preproc.warp.cleanup     = 1;
            matlabbatch{gi}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
            matlabbatch{gi}.spm.spatial.preproc.warp.affreg      = 'mni';
            matlabbatch{gi}.spm.spatial.preproc.warp.fwhm        = 0;
            matlabbatch{gi}.spm.spatial.preproc.warp.samp        = 3;
            matlabbatch{gi}.spm.spatial.preproc.warp.write       = [0 0];
            gi = gi + 1;
        end
        
        %-------------------------------
        %Do skull strip
        if do_skull
            Vfnames      = strvcat(struc_file,c1_file,c2_file);
            matlabbatch{gi}.spm.util.imcalc.input            = cellstr(Vfnames);
            matlabbatch{gi}.spm.util.imcalc.output           = skullstrip_name;
            matlabbatch{gi}.spm.util.imcalc.outdir           = {st_dir};
            matlabbatch{gi}.spm.util.imcalc.expression       = 'i1.*((i2+i3)>0.2)';
            matlabbatch{gi}.spm.util.imcalc.options.dmtx     = 0;
            matlabbatch{gi}.spm.util.imcalc.options.mask     = 0;
            matlabbatch{gi}.spm.util.imcalc.options.interp   = 1;
            matlabbatch{gi}.spm.util.imcalc.options.dtype    = 4;
            gi = gi + 1;
        end
        %-------------------------------
        %Dartel norm to template
        if do_norm
            matlabbatch{gi}.spm.tools.dartel.warp1.images = {cellstr(rc1_file),cellstr(rc2_file)};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.rform = 0;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(1).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(1).K = 0;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(1).template = {[template_path 'Template_1_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(2).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(2).K = 0;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(2).template = {[template_path 'Template_2_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(3).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(3).K = 1;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(3).template = {[template_path 'Template_3_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(4).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(4).K = 2;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(4).template = {[template_path 'Template_4_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(5).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(5).K = 4;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(5).template = {[template_path 'Template_5_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(6).its = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(6).K = 6;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.param(6).template = {[template_path 'Template_6_IXI555_MNI152.nii']};
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
            matlabbatch{gi}.spm.tools.dartel.warp1.settings.optim.its = 3;
            gi = gi + 1;
        end
        %-------------------------------
        %Get backwards deformations
        if do_back
            matlabbatch{gi}.spm.util.defs.comp{1}.dartel.flowfield = {u_rc1_file};
            matlabbatch{gi}.spm.util.defs.comp{1}.dartel.times     = [1 0];
            matlabbatch{gi}.spm.util.defs.comp{1}.dartel.K         = 6;
            matlabbatch{gi}.spm.util.defs.comp{1}.dartel.template  = {''};
            matlabbatch{gi}.spm.util.defs.out{1}.savedef.ofname    = 'backwards';
            matlabbatch{gi}.spm.util.defs.out{1}.savedef.savedir.saveusr = {st_dir};
            gi = gi + 1;
        end
        
        %-------------------------------
        %Create warped T1 and mean EPI
        if do_warp
            matlabbatch{gi}.spm.tools.dartel.crt_warped.flowfields = cellstr(strvcat(u_rc1_file,u_rc1_file,u_rc1_file,u_rc1_file));
            matlabbatch{gi}.spm.tools.dartel.crt_warped.images = {cellstr(strvcat(mean_file,strip_file,c1_file,c2_file))};
            matlabbatch{gi}.spm.tools.dartel.crt_warped.jactransf = 0;
            matlabbatch{gi}.spm.tools.dartel.crt_warped.K = 6;
            matlabbatch{gi}.spm.tools.dartel.crt_warped.interp = 1;
            gi = gi + 1;
        end
        %-------------------------------
        %Create smoothed skullstrip
        if do_sm_skull
            skern = 3;
            matlabbatch{gi}.spm.spatial.smooth.data   = cellstr(strip_file);
            matlabbatch{gi}.spm.spatial.smooth.fwhm   = repmat(skern,1,3);
            matlabbatch{gi}.spm.spatial.smooth.prefix = ['s' num2str(skern)];
            gi = gi + 1;
        end
    end
    if ~isempty(matlabbatch)
        run_matlab(np, matlabbatch, check);
    end
end

if do_avg_norm
    matlabbatch = [];
    all_wskull_files = [];
    all_wmean_files  = [];
    all_wc1_files    = [];
    
    for g = 1:size(all_subs,2)
        name = sprintf('Sub%02.2d',all_subs(g));
        strip_file        = [base_dir name filesep 'T1' filesep skullstrip_name];
        wskull_file       = ins_letter(strip_file,'w');
        
        mean_file        = [base_dir name filesep 'T1' filesep mean_func_name];
        wmean_file       = ins_letter(mean_file,'w');


        
        st_dir       = [base_dir name filesep 'T1' filesep];
        struc_file   = spm_select('FPList', st_dir, struc_templ);
        wc1_file     = ins_letter(struc_file,'wc1');

        all_wskull_files  = strvcat(all_wskull_files,wskull_file);
        all_wmean_files   = strvcat(all_wmean_files,wmean_file);
        all_wc1_files     = strvcat(all_wc1_files,wc1_file);
    end
    
    matlabbatch{1}.spm.util.imcalc.input = cellstr(all_wskull_files);
    matlabbatch{1}.spm.util.imcalc.output = 'mean_wskull';
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(base_dir);
    matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    matlabbatch{2} = matlabbatch{1};
    matlabbatch{2}.spm.util.imcalc.input = cellstr(all_wmean_files);
    matlabbatch{2}.spm.util.imcalc.output = 'mean_wmean';
    
    matlabbatch{3} = matlabbatch{1};
    matlabbatch{3}.spm.util.imcalc.input = cellstr(all_wc1_files);
    matlabbatch{3}.spm.util.imcalc.output = 'mean_wc1';

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end

function chuckCell = splitvect(v, n)
% Splits a vector into number of n chunks of  the same size (if possible).
% In not possible the chunks are almost of equal size.
%
% based on http://code.activestate.com/recipes/425044/
chuckCell  = {};
vectLength = numel(v);
splitsize  = 1/n*vectLength;
for i = 1:n
    idxs = [floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1]+1;
    chuckCell{end + 1} = v(idxs);
end


function f_files = create_func_files(s_dir,f_templ,n_files)
for i=1:n_files
    f_files{i} = [s_dir filesep f_templ ',' num2str(i)];
end


function out = ins_letter(pscan,letter_start,letter_end)
if nargin <3
    letter_end = [];
end
for a=1:size(pscan,1)
    [p , f, e] = fileparts(pscan(a,:));
    out(a,:) = [p filesep letter_start f letter_end e];
end

function run_matlab(np, matlabbatch, check)

spm_path          = fileparts(which('spm')); %get spm path
mat_name          = which(mfilename);
[~,mat_name,~]    = fileparts(mat_name);


fname = [num2str(np) '_' mat_name '.mat'];

save([num2str(np) '_' mat_name],'matlabbatch');
lo_cmd  = ['clear matlabbatch;load(''' fname ''');'];
ex_cmd  = ['addpath(''' spm_path ''');spm(''defaults'',''FMRI'');spm_jobman(''initcfg'');spm_jobman(''run'',matlabbatch);'];
end_cmd = ['delete(''' fname ''');'];
if ~check
    system(['start matlab.exe -nodesktop -nosplash  -logfile ' num2str(np) '_' mat_name '.log -r "' lo_cmd ex_cmd end_cmd 'exit"']);
end

function out = dis_fname(pfname)

for i = 1:size(pfname,1)
    [p fname ext] = fileparts(pfname(i,:));
    
    
    us = strfind(fname,'_');
    hy = strfind(fname,'-');
    
    out(i).path   = p;
    out(i).pre    = fname(1:us-1);
    out(i).sub    = str2num(fname(us+1:hy(1)-1));
    out(i).series = str2num(fname(hy(1)+1:hy(2)-1));
    out(i).ind1   = str2num(fname(hy(2)+1:hy(3)-1));
    out(i).ind2   = str2num(fname(hy(3)+1:hy(4)-1));
    out(i).count  = str2num(fname(hy(4)+1:end));
    out(i).ext    = ext;
end
