function [constant,crossmod] = analysis_crossmod_preproc_beam(constant)

curdir = cd;

%% Cross-sectional Modeller

cd('cross_mod')
[crossmod] = cross_mod_beam(constant);
cd(curdir)
fprintf(' Cross-sectional modeller ;')
