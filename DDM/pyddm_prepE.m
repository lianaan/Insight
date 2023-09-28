clear all; close all;
load('alldata_E2.mat')

Nsubj = size(alldata,1);
Ntrials = length(alldata(1,1).stims); % 121
Ncond = size(alldata,2);


delimiter = ',';

for si = 1: 1: Nsubj
    for ci = 1: Ncond
    
    file_name_i = ['mae', '_s_', num2str(si),'_c_', num2str(ci), '.csv'];
    
    fid = fopen(file_name_i,'a');   

    fprintf(fid,' %s','participant'); 
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','rt');  
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','coh');  
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','correct');  
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','trgchoice');  % left is 0 and right is 1
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','left_is_correct');  
    fprintf(fid, [delimiter]);
    
    fprintf(fid,' %s','highreward');  % if the participant chose right (1)
    fprintf(fid, [delimiter]);
    
    fprintf(fid, '\n'); % move to the next line
    
    for ti = 1: Ntrials
        fprintf(fid, '%i',1);
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '%0.4f',alldata(si,ci).resp_times(ti));
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '%0.4f',abs(alldata(si,ci).stims(ti)));
        fprintf(fid, [delimiter]);
        
        correct_val = [(alldata(si,ci).stims(ti)>=0 & alldata(si,ci).resp(ti) == 1) |(alldata(si,ci).stims(ti)<=0 & alldata(si,ci).resp(ti) == 0)];  
        fprintf(fid, '%i',correct_val);
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '%i',alldata(si,ci).resp(ti)+1);
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '%i',alldata(si,ci).stims(ti)<=0);
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '%i',alldata(si,ci).resp(ti)==1);
        fprintf(fid, [delimiter]);
        
        fprintf(fid, '\n');
        
    end
    fclose(fid);
    end
end



