function []=loadStructures(n,step,type)
% idx=[11 50 12 51 17 53 18 54];
idx=[17 53 18 54];
if nargin<1
    ls '/extra/hermes/scratch/dbreton/ADNI/freesurfer_5.1/results'
else
    addpath '/extra/disco1/miaTools/BrainTools';
    if nargin<3, type='original'; end;
    
    FS=read_FS(['/extra/hermes/scratch/dbreton/ADNI/freesurfer_5.1/results/RID_NOR_' num2str(n,'%04d') '/' ],type);
    if nargin<2, step=1; end;
    for i=1:size(idx,1)
            disp(['STRUCTURE ' num2str(i)]);
            IS=isosurface (FS.fill (FS{idx(i)}),0.5+1e-5);
            Miso{i}.xyz=IS.vertices;
            Miso{i}.tri=IS.faces;
            M{i}=AFRONT_isosurface (FS.fill (FS{idx(i)}), 0.5+1e-5 ,['-min_step ' num2str(step)]);

   end;
   save(['ST' num2str(n,'%04d')],'M','Miso');
    
end;
