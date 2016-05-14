clc
clear all
close all

dir = pwd;

fn = fullfile([dir, '\compare']);
if ~ exist(fn, 'dir')
   mkdir(fn)
end

copyfile(strcat(dir, '\FAST_WInDS_rigid_noDS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_WInDS_rigid_noDS.outb'));

copyfile(strcat(dir, '\FAST_WInDS_rigid_DS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_WInDS_rigid_DS.outb'));

copyfile(strcat(dir, '\FAST_WInDS_elastic_noDS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_WInDS_elastic_noDS.outb'));

copyfile(strcat(dir, '\FAST_WInDS_elastic_DS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_WInDS_elastic_DS.outb'));



copyfile(strcat(dir, '\FAST_BEM_rigid_noDS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_BEM_rigid_noDS.outb'));

copyfile(strcat(dir, '\FAST_BEM_rigid_DS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_BEM_rigid_DS.outb'));

copyfile(strcat(dir, '\FAST_BEM_elastic_noDS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_BEM_elastic_noDS.outb'));

copyfile(strcat(dir, '\FAST_BEM_elastic_DS\Test24.outb'), ...
         strcat(dir, '\compare\Test24_FAST_BEM_elastic_DS.outb'));


disp('Copied...')     


