clear

if verLessThan('matlab','8.2')
    matlabpool open
else
    poolobj = parpool;
end;

% Array containing the weights
wArr = zeros(32,1);

% Uniform
% for sizes=[60 64]
%     % Init
%     arr = AntArray(zeros(sizes), 60500, [], 0.84);
%     arr = arr.setName(['fullx_' mat2str(sizes)]);
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
%     arr = arr.setComments(sprintf('Elements spacing: 0.84$\\lambda$'));
%     
%     % Plots
%     arr.plotAntArray();
%     arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
%     arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
%     arr = arr.genPattern(11000, 3000, 'XY', 30);
%     arr = arr.genPattern([], [], 'XY-BW');
% %     if sizes == 60
% %         wArr(1,1) = arr.weight('XY');
% %         wArr(1,2) = arr.weight('theta', pi/4);
% %     else
% %         wArr(2,1) = arr.weight('XY');
% %         wArr(2,2) = arr.weight('theta', pi/4);
% %     end;
%     arr = arr.E_strength(15000, 0, 0, 500);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% end;

% Squares
% for rem_els=0:15
%     disp(['Squares, iteration ' num2str(rem_els) ' of 15']);
%     % Init
%     arr = AntArray(zeros(64), 60500, [], 0.84);
%     arr = arr.setName(['sq_' mat2str(rem_els)]);
%     arr = arr.setMax('XY', 30);
%     arr = arr.setMax('YZ', 30);
%     arr = arr.setMax('E', 25);
%     arr = arr.setMin('XY', -60);
%     arr = arr.setMin('YZ', -60);
%     arr = arr.setMin('E', -15);
% 
%     % Create elements' pattern
%     tmp = zeros(64);
%     for i=16:32-rem_els
%         tmp(i, i:end-i+1)=1;
%     end;
% 
%     tmp(tmp(end:-1:1,:)==1)=1;
%     tmp(tmp'==1)=1;
% 
%     for j=1:16-rem_els
%         for i=1:32
%             ln = i+2*(j-1);
%             rw = 32-i+1;
%             if ln <= 32 && rw <=32
%                 tmp(ln,rw)=1;
%             end;
%         end;
%     end;
% 
%     tmp(tmp'==1)=1;
%     tmp(tmp(end:-1:1,:)==1)=1;
%     tmp(tmp(:,end:-1:1)==1)=1;
% 
%     arr = arr.adaptArray(tmp, 90000, 0, 0);
%     
%     el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
%     
%     arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
%         'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));
% 
%     % Plots
%     arr.plotAntArray();
%     arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
%     arr = arr.genPattern([], [], 'theta-BW', [], 0);
%     arr = arr.genPattern(11000, 3000, 'XY', 30);
%     arr = arr.genPattern([], [], 'XY-BW');
%     arr = arr.E_strength(15000, 0, 0, 500);
%     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
%     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% end;

% for rem_els=0:15
%     disp(['Squares 2, iteration ' num2str(rem_els) ' of 15']);
%     % Init
%     arr = AntArray(zeros(64), 60500, [], 0.84);
%     arr = arr.setName(['sq2_' mat2str(rem_els)]);
%     arr = arr.setMax('XY', 30);
%     arr = arr.setMax('YZ', 30);
%     arr = arr.setMax('E', 25);
%     arr = arr.setMin('XY', -60);
%     arr = arr.setMin('YZ', -60);
%     arr = arr.setMin('E', -15);
% 
%     % Create elements' pattern
%     tmp = zeros(64);
%     for i=16:31-rem_els
%         tmp(i, i:end-i+1)=1;
%     end;
% 
%     tmp(tmp(end:-1:1,:)==1)=1;
%     tmp(tmp'==1)=1;
% 
%     for j=1:16-rem_els
%         for i=1:32
%             ln = i+2*(j-1);
%             rw = 32-i+1;
%             if ln <= 32 && rw <=32
%                 tmp(ln,rw)=1;
%             end;
%         end;
%     end;
% 
%     tmp(tmp'==1)=1;
%     tmp(tmp(end:-1:1,:)==1)=1;
%     tmp(tmp(:,end:-1:1)==1)=1;
% 
%     arr = arr.adaptArray(tmp, 90000, 0, 0);
%     
%     el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
%     
%     arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
%         'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));
% 
%     % Plots
% %     arr.plotAntArray();
% %     arr = arr.genPattern(11000, 3000, 'theta', 30, pi/4);
% %     arr = arr.genPattern([], [], 'theta-BW', [], pi/4);
%     wArr(rem_els+1,3) = arr.weight('theta', pi/4);
% %     arr = arr.genPattern(11000, 3000, 'XY', 30);
% %     arr = arr.genPattern([], [], 'XY-BW');
%     wArr(rem_els+1,4) = arr.weight('XY');
% %     arr = arr.E_strength(15000, 0, 0, 500);
% %     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
% %     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% end;
% 
% % Triangles
for rem_els=9
    disp(['Triangles, iteration ' num2str(rem_els) ' of 15']);
    % Init
    arr = AntArray(zeros(60), 60500, [], 0.84);
    arr = arr.setName(['tr_' num2str(rem_els)]);
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
    
%     arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
%         'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));

    % Plots
    arr.plotAntArray();
% %     arr = arr.genPattern(11000, 3000, 'theta', 30, 0);
% %     arr = arr.genPattern([], [], 'theta-BW', [], 0);
%     wArr(rem_els+1,5) = arr.weight('theta', 0);
% %     arr = arr.genPattern(11000, 3000, 'XY', 30);
% %     arr = arr.genPattern([], [], 'XY-BW');
%     wArr(rem_els+1,6) = arr.weight('XY');
% %     arr = arr.E_strength(15000, 0, 0, 500);
% %     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
% %     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% end;

% Triangles 2
% for rem_els=0:15
%     disp(['Triangles 2, iteration ' num2str(rem_els) ' of 15']);
%     % Init
    arr = AntArray(zeros(60), 60500, [], 0.84);
    arr = arr.setName(['tr2_' num2str(rem_els)]);
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

    tmp(tmp'==1)=1;
    tmp(tmp(end:-1:1,:)==1)=1;
    tmp(tmp(:,end:-1:1)==1)=1;

    arr = arr.adaptArray(tmp, 90000, 0, 0);
    
    el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
%     
%     arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
%         'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));
% 
%     % Plots
    arr.plotAntArray();
% %     arr = arr.genPattern(11000, 3000, 'XY', 30);
% %     arr = arr.genPattern([], [], 'XY-BW');
%     wArr(rem_els+1,7) = arr.weight('XY');
% %     arr = arr.E_strength(15000, 0, 0, 500);
% %     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
% %     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
end;

% Circles
% for rem_els=0:31
%     disp(['Circles, iteration ' num2str(rem_els) ' of 31']);
%     % Init
%     arr = AntArray(zeros(64), 60500, [], 0.84);
%     arr = arr.setName(['circ_' num2str(rem_els)]);
%     arr = arr.setMax('XY', 30);
%     arr = arr.setMax('YZ', 30);
%     arr = arr.setMax('E', 25);
%     arr = arr.setMin('XY', -60);
%     arr = arr.setMin('YZ', -60);
%     arr = arr.setMin('E', -15);
% 
%     % Create elements' pattern
%     tmp = drawCircle(64, 32, rem_els);
% 
%     arr = arr.adaptArray(tmp, 90000, 0, 0);
%     
%     el_ratio = num2str(length(find(tmp~=0))/numel(tmp)*100,3);
%     arr = arr.setComments(sprintf([num2str(rem_els) ' lines removed\n' ...
%         'Elements spacing: 0.84$\\lambda$\n\\# of elements: ' el_ratio '\\%%']));
% 
%     % Plots
%     arr.plotAntArray();
% %     arr = arr.genPattern(11000, 3000, 'XY', 30);
% %     arr = arr.genPattern([], [], 'XY-BW');
%     wArr(rem_els+1,1) = arr.weight('XY');
% %     arr = arr.E_strength(15000, 0, 0, 500);
% % %     arr = arr.genPattern([1 2:2:10]*1000, 3000, 'YZ', 30);
% % %     arr = arr.genPattern([1 2:2:10]*1000, [], 'YZ-BW');
% %     arr = arr.genPattern(10*1000, 3000, 'YZ', 30);
% %     arr = arr.genPattern(10*1000, [], 'YZ-BW');
% 
% end;
% 
% export_dat(wArr, 'weights_circ');

if verLessThan('matlab','8.2')
    matlabpool close
else
   delete(poolobj);
end;