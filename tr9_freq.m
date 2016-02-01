clear

if verLessThan('matlab','8.2')
    matlabpool open
else
    poolobj = parpool;
end;


% Uniform
% Init
% sizes = 60;
% space = 0.84;
% for freq=61.2:0.2:61.6
%     arr = AntArray(zeros(sizes), freq*1e3, [], space);
%     arr = arr.setName(['fullx' mat2str(sizes) '_f' mat2str(freq*10)]);
%     arr = arr.setNormFreq(60500);
%     arr = arr.setMax('XY', 30);
%     arr = arr.setMax('YZ', 30);
%     arr = arr.setMax('E', 25);
%     arr = arr.setMin('XY', -60);
%     arr = arr.setMin('YZ', -60);
%     arr = arr.setMin('E', -15);
% 
%     % Create elements' pattern
%     arr = arr.adaptArray(ones(sizes), 90000, 0, 0);
% 
%     arr = arr.setComments(sprintf(['Elements spacing: ' mat2str(space) ...
%         '$\\lambda$\nFrequency: ' mat2str(freq) 'GHz']));
% 
%     % Plots
%     arr.plotAntArray();
%     arr = arr.genPattern(11000, 3000, 'XY', 30);
%     arr = arr.genPattern([], [], 'XY-BW');
% %     arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
% %     arr = arr.genPattern([], [], 'theta-BW', [], 0);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% end;

% Element spacing influence
for freq=61.6:0.2:61.6
    space = 0.84;
    rem_els = 9;
    % Init
    arr = AntArray(zeros(60), freq*1e3, [], space);
    arr = arr.setName(['tr' mat2str(rem_els) '_f' mat2str(freq*10)]);
    arr = arr.setNormFreq(60500);
    arr = arr.setMax('XY', 30);
    arr = arr.setMax('YZ', 30);
    arr = arr.setMax('E', 25);
    arr = arr.setMin('XY', -60);
    arr = arr.setMin('YZ', -60);
    arr = arr.setMin('E', -15);

    % Create elements' pattern
    tmp = zeros(60);
    len = size(arr.M,1)/4;
    d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
    k_max = 30 - round(rem_els/d_ratio);
    k_lim = round(k_max*d_ratio);
    for k=1:k_max
        i_max = size(arr.M,1)-len+2-k;
        if k <= k_lim
            tmp(len-1+k, 3+round(k/tan(pi/3)):end-2-round(k/tan(pi/3)))=1;
        end;
        for i=min(k,max(k_lim - 1,1)):i_max
           tmp(len-1+i, 3+round((i+k-1)/tan(pi/3)))=1; 
        end
    end;

    tmp(tmp(end:-1:1,:)==1)=1;
    tmp(tmp(:,end:-1:1)==1)=1;

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
    
    arr = arr.setComments(sprintf(['Elements spacing: ' mat2str(space) ...
        '$\\lambda$\nFrequency: ' mat2str(freq) 'GHz']));

    % Plots
    arr.plotAntArray();
    arr = arr.genPattern(11000, 3000, 'XY', 30);
    arr = arr.genPattern([], [], 'XY-BW');
    arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
    arr = arr.genPattern([], [], 'theta-BW', [], 0);
    arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
    arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

if verLessThan('matlab','8.2')
    matlabpool close
else
   delete(poolobj);
end;