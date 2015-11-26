function [] = MTF_export_figures(final_tracks, tracks_tail, seq_name)
% MTF_export_figures(final_tracks, tail_tracks)
%
% Plots the tracks resulting from the MTF tracker. The default plots are:
%
% X vs time/frame for 4 paws.
% X vs Z for FR FL
% X vs Z for HR HL
% Z vs time for FR FL
% Z vs time for HR HL
% Y vs X for tail
% Z vs X for tail

[~,seq_name,~] = fileparts(seq_name);

% plotting tracks:
Colors = ['r';'m';'b';'c';'g'];
Legend = {'FR','HR','FL','HL','T'};
N_images = size(final_tracks,3);

% X vs time/frame for 4 paws.
figure()
hold on
for i_paw = 1:4
    plot(1:N_images,squeeze(final_tracks(1,i_paw,:)),'-','Color',Colors(i_paw,:),'LineWidth',2);
end
legend(Legend(1:4),'Location','NorthWest');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
saveas(gcf,sprintf('%s_x_vs_t.png',seq_name));
close(gcf)

% Z vs time for FR FL
% Z vs time for HR HL
pairs = [1 2;3 4];
figure()
for i_pairs = 1:2
    subplot(2,1,i_pairs)
    hold on
    plot(1:N_images,squeeze(final_tracks(3,pairs(i_pairs,1),:)),'-','Color',Colors(pairs(i_pairs,1),:),'LineWidth',2);
    plot(1:N_images,squeeze(final_tracks(3,pairs(i_pairs,2),:)),'-','Color',Colors(pairs(i_pairs,2),:),'LineWidth',2);
    legend(Legend(pairs(i_pairs,:)),'Location','NorthWest');
    axis ij tight
end
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
saveas(gcf,sprintf('%s_Time_vs_Z.png',seq_name));
close(gcf)

% Y vs X for tail
type = {'X vs Z','X vs Y'};
figure()
idx = {[1 2],[1 3]};
for i_tail = 1:2
    subplot(2,1,i_tail)
    plot(squeeze(tracks_tail(idx{i_tail}(1),:,:))',squeeze(tracks_tail(idx{i_tail}(2),:,:))','-','LineWidth',2);
    legend(type{i_tail},'Location','NorthEast');
%     axis ij equal tight  
end
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
saveas(gcf,sprintf('%s_%s.png',seq_name,Legend{5}));
drawnow;
close(gcf)
