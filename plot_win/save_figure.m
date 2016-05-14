clc
clear
close all
 
%% User inputs
wind = '_above';   
colunms={'RotSpeed','RotThrust', 'RotTorq','GenPwr','Wave1Elev','PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw'};
test_num = '24';
output_dir = '/figures'; % under 'compare' 

%.................................................................................................................
% Load files
    
dir = pwd;   
full_dir = strcat(dir,'/compare');
    
[Channels8, ChanName8, ChanUnit8, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_WInDS_rigid_noDS.outb'));
[Channels7, ChanName7, ChanUnit7, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_WInDS_rigid_DS.outb'));
[Channels6, ChanName6, ChanUnit6, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_BEM_rigid_noDS.outb'));
[Channels5, ChanName5, ChanUnit5, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_BEM_rigid_DS.outb'));   

[Channels4, ChanName4, ChanUnit4, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_WInDS_elastic_noDS.outb'));
[Channels3, ChanName3, ChanUnit3, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_WInDS_elastic_DS.outb'));
[Channels2, ChanName2, ChanUnit2, FileID] = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_BEM_elastic_noDS.outb'));
[Channels,  ChanName,  ChanUnit, FileID]  = ReadFASTbinary(strcat(full_dir, '/Test',test_num,'_FAST_BEM_elastic_DS.outb'));


fn = fullfile([full_dir, output_dir]);
if ~ exist(fn, 'dir')
   mkdir(fn)
end


%% Plot

for i=1:length(colunms)
    figure(i) 
        for j=1:length(ChanName)
            if strcmp(ChanName(j),colunms(i)) 
               index = j; 
            end
        end    
    
    if strcmp(colunms(i), 'Wave1Elev')
        plot(Channels(:,1)', Channels(:,index)','b')
        
    else     
        plot(Channels(:,1)', Channels(:,index)','r')
        hold on

        plot(Channels2(:,1)', Channels2(:,index)','k');
        plot(Channels3(:,1)', Channels3(:,index)','g');
        plot(Channels4(:,1)', Channels4(:,index)','b');

        plot(Channels5(:,1)', Channels5(:,index)','r-.')
        plot(Channels6(:,1)', Channels6(:,index)','k-.');
        plot(Channels7(:,1)', Channels7(:,index)','g-.');
        plot(Channels8(:,1)', Channels8(:,index)','b-.');    

        hleg1 = legend('BEM elastic DS','BEM elastic no DS','WInDS elastic DS','WInDS elastic no DS', 'BEM rigid DS','BEM rigid no DS','WInDS rigid DS','WInDS rigid no DS');
        set(hleg1,'Location','NorthEast')
    end
    
    xlim([0 max(Channels(:,1))])
    grid on

    xlabel('Time (sec)')
    ylabel([colunms(i),ChanUnit(index) ])
    % title(['Rotor Speed of Test ' TestNum], 'FontSize', 12)

    name = colunms(i);
    name = name{:};
    saveas(figure(i), [full_dir, output_dir, '/Test',test_num,'_all', wind,'_',name,'.fig'])
    saveas(figure(i), [full_dir, output_dir, '/Test',test_num,'_all', wind,'_',name],'epsc')    
 end
