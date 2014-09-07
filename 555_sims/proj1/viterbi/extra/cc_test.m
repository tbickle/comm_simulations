close all; clear all; clc;
load cc_test.mat

K = 3;
G = [1 1 1;1 1 0;1 1 1];

ns = [0 0]; % original state
o_data = [1 1 0 0 0]
fprintf('t=%d  out=n/a  state=%d%d\n',0,ns(1,1),ns(1,2)); % test
for i=1:length(o_data)*0+5
    % output while transitioning to next state
    o_out = mod([o_data(i) ns]*G,2);
    if(i==1) signal=o_out;
    else signal = horzcat(signal,o_out);
    end
    % next state
    ns = [o_data(i) ns(1)];
    % display
    fprintf('t=%d  out=%d%d%d  state=%d%d\n',i,o_out(1,1),o_out(1,2),o_out(1,3),ns(1,1),ns(1,2)); % test
end
signal

% MODULATE
% CHANNEL
% DEMODULATE

state = [0 0]; % assumed initial state
% bad_data = [1 0 1 0 0 1 1 0 1 1 1 1 0 0 0];
data = signal;
b_est = zeros(1,length(data)/K);
fprintf('t=%d state=%d%d\n',0,state(1,1),state(1,2)); % test
for i=1:length(data)/K
    rx = data(1,K*(i-1)+1:K*(i-1)+K); % extract data from rx'd signal
    td = [horzcat(0,state);horzcat(1,state)];
    [value index] = min(sum(abs(mod(td*G,2)-[rx;rx]),2));
    state = td(index,1:2);
    b_est(1,i) = state(1,1);
    fprintf('t=%d state=%d%d\n',i,state(1,1),state(1,2)); % test
end
b_est

error = sign(abs(o_data-b_est))
