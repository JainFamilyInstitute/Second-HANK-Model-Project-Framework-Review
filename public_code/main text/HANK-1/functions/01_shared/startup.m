%
% Get IP-address from DOS-script
% [cmd_status,ipv4addr]=dos('ipaddr.cmd');
%
% Configure IP-address for Parallel Computing Toolbox client session
pctconfig('hostname','131.220.108.212');
% pctconfig('hostname',ipv4addr);

%
% Add all other commands
% userpath('...');