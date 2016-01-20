clear

if verLessThan('matlab','8.2')
    matlabpool open
else
    poolobj = parpool;
end;

weightArr = zeros(9,2);

% % Uniform
% % Init
% sizes = 60;
% arr = AntArray(zeros(sizes), 60500, [], 0.84);
% arr = arr.setName(['fullx' mat2str(sizes)]);
% arr = arr.setMax('XY', 30);
% arr = arr.setMax('YZ', 30);
% arr = arr.setMax('E', 25);
% arr = arr.setMin('XY', -60);
% arr = arr.setMin('YZ', -60);
% arr = arr.setMin('E', -15);
% 
% % Create elements' pattern
% arr = arr.adaptArray(ones(sizes), 90000, 0, 0);
% 
% arr = arr.setComments(sprintf('Elements spacing: 0.84$\\lambda$'));
% 
% % Plots
% arr.plotAntArray();
% for ang=1:2:11
%     disp(['Current angle: ' rats(ang/24) 'pi']);
%     arr = arr.genPattern(11000, 3000, 'theta', 30, ang*pi/24);
%     arr = arr.genPattern([], [], 'theta-BW', [], ang*pi/24);
% end;


% % Triangles
% for rem_els=9
%     % Init
%     arr = AntArray(zeros(60), 60500, [], 0.84);
%     arr = arr.setName(['tr' num2str(rem_els)]);
%     arr = arr.setMax('XY', 30);
%     arr = arr.setMax('YZ', 30);
%     arr = arr.setMax('E', 25);
%     arr = arr.setMin('XY', -60);
%     arr = arr.setMin('YZ', -60);
%     arr = arr.setMin('E', -15);
% 
%     % Create elements' pattern
%     tmp = zeros(60);
%     len = size(arr.M,1)/4;
%     d_ratio = abs(1-tan(pi/3))/sqrt((tan(pi/3))^2-1);
%     k_max = 30 - round(rem_els/d_ratio);
%     k_lim = round(k_max*d_ratio);
%     for k=1:k_max
%         i_max = size(arr.M,1)-len+2-k;
%         if k <= k_lim
%             tmp(len-1+k, 3+round(k/tan(pi/3)):end-2-round(k/tan(pi/3)))=1;
%         end;
%         for i=min(k,max(k_lim - 1,1)):i_max
%            tmp(len-1+i, 3+round((i+k-1)/tan(pi/3)))=1; 
%         end
%     end;
% 
%     tmp(tmp(end:-1:1,:)==1)=1;
%     tmp(tmp(:,end:-1:1)==1)=1;
% 
%     arr = arr.adaptArray(tmp, 90000, 0, 0);
%     
%     el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
%     
%     arr = arr.setComments(sprintf('Elements spacing: 0.84$\\lambda$'));
% 
%     % Plots
%     arr.plotAntArray();
%     for ang=1:2:11
%         disp(['Current angle: ' rats(ang/24) 'pi']);
%         arr = arr.genPattern(11000, 3000, 'theta', 30, ang*pi/24);
%         arr = arr.genPattern([], [], 'theta-BW', [], ang*pi/24);
%     end;
% end;

% Element spacing influence
for space=0.84:0.02:1
    % Init
    arr = AntArray(zeros(60), 60500, [], space);
    arr = arr.setName(['tr8_' mat2str(space*100)]);
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
    k_max = 30 - round(8/d_ratio);
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
        '$\\lambda$']));

    % Plots
    iter = round((space - 0.84)/0.02 + 1);
    
    arr.plotAntArray();
%     arr = arr.genPattern(11000, 3000, 'XY', 30);
%     arr = arr.genPattern([], [], 'XY-BW');
    weightArr(iter, 1) = arr.weight('XY');
%     arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
%     arr = arr.genPattern([], [], 'theta-BW', [], 0);
    weightArr(iter, 2) = arr.weight('theta', 0);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

for space=0.84:0.02:1
    % Init
    arr = AntArray(zeros(60), 60500, [], space);
    arr = arr.setName(['tr9_' mat2str(space*100)]);
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
    k_max = 30 - round(9/d_ratio);
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
        '$\\lambda$']));

    % Plots
    iter = round((space - 0.84)/0.02 + 1);
    
    arr.plotAntArray();
%     arr = arr.genPattern(11000, 3000, 'XY', 30);
%     arr = arr.genPattern([], [], 'XY-BW');
    weightArr(iter, 3) = arr.weight('XY');
%     arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
%     arr = arr.genPattern([], [], 'theta-BW', [], 0);
    weightArr(iter, 4) = arr.weight('theta', 0);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;


export_dat(weightArr, 'weights');

if verLessThan('matlab','8.2')
    matlabpool close
else
   delete(poolobj);
end;