function [T5, T6] = soft_sync_cam_nlx(event_file)
% IEEE-1588 clock synchronization
% October 27, 2020
% Author Shahin G. Lashkari

addpath('../jumping');
[Timestamps,~,~,EventIDs,TTLs] = readevent(event_file);

T5 = Timestamps(TTLs==5 & EventIDs==11);
T6 = Timestamps(TTLs==6 & EventIDs==11);

end