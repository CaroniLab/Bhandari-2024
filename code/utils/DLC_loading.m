%% DLC node
dlc_file = dir(fullfile(path, "*DLC*filtered.csv")).name;
opts = detectImportOptions((fullfile(path, dlc_file)));
opts.DataLines = 4;
opts.VariableNamesLine = 3;
opts.VariableDescriptionsLine= 2;
coords = readtable((fullfile(path, dlc_file)), opts);
colnames= ["ts"];
bodyparts = ["Snout", "tailbase", "bodycentre", "object1", "object2"];
for ii = 1:length(bodyparts)
    colnames = [colnames, strcat(bodyparts(ii), "_x"), strcat(bodyparts(ii), "_y"), strcat(bodyparts(ii), "_likelihood")];
end
coords.Properties.VariableNames = colnames;
coords.ts = camera_ts;
%%
%FOR_stop = FOR_start + 3000;
coords_FOR = coords((coords.ts >= FOR_start) & (coords.ts <= FOR_stop), :); 
camera_ts_FOR = camera_ts((camera_ts>= FOR_start) & (camera_ts<= FOR_stop));
figure
plot(coords_FOR.Snout_x, coords_FOR.Snout_y)
title("Before cleanup of data")
%% Medfilt and clean up tracking data
if manual_obj == 1
    objects_loc = csvread(fullfile(path, 'manual_object_location.csv'));
    coords_FOR.("object1_x")(:) =  objects_loc(1,1);
    coords_FOR.("object1_y")(:) =  objects_loc(1,2);
    coords_FOR.("object2_x")(:) =  objects_loc(2,1);
    coords_FOR.("object2_y")(:) =  objects_loc(2,2);   
end

%% 
figure
p_thresh = 0.9;
for ii = 1:length(bodyparts)
    bodypart = bodyparts(ii);
    lh = coords_FOR.(bodypart+"_likelihood");
    x = coords_FOR.(bodypart+"_x");
    y = coords_FOR.(bodypart+"_y");
    
    x(lh < p_thresh) = NaN;
    y(lh < p_thresh) = NaN;
    coords_FOR.(bodypart+"_x") = fillmissing(x, 'linear');
    coords_FOR.(bodypart+"_y") = fillmissing(y, 'linear');
    coords_FOR.(bodypart+"_x") = medfilt1(coords_FOR.(bodypart+"_x"), 15);
    coords_FOR.(bodypart+"_y") = medfilt1(coords_FOR.(bodypart+"_y"), 15);
    

    
    plot(coords_FOR.(bodypart+"_x"), coords_FOR.(bodypart+"_y"))
    hold on
end
legend(bodyparts)
scatter(coords_FOR.("object1_x")(50),coords_FOR.("object1_y")(50))
scatter(coords_FOR.("object2_x")(50),coords_FOR.("object2_y")(50))
title("After cleanup of data")  
saveas(gcf, fullfile(path, 'plots_DLC', sprintf('animal_trajectory.png')))
%% Distance between Snout and object
snout_obj1 =vecnorm(coords_FOR{:, ["Snout_x", "Snout_y"]} - coords_FOR{:, ["object1_x", "object1_y"]}, 2, 2);
snout_obj2 =vecnorm(coords_FOR{:, ["Snout_x", "Snout_y"]} - coords_FOR{:, ["object2_x", "object2_y"]}, 2, 2);
figure
plot(camera_ts_FOR, snout_obj1, 'r')
hold on
yyaxis right
plot(camera_ts_FOR, snout_obj2, 'k')
legend("d from obj 1", "d from obj 2")

%%
%uses thresholds set in analysis script

%figure
[dist_approach_obj2, obj2_app_DLC_frames] = findpeaks(-snout_obj2, 'MinPeakHeight', -abs(thresh), 'MinPeakDistance',200, 'MinPeakProminence',15);
[dist_approach_obj1, obj1_app_DLC_frames] = findpeaks(-snout_obj1, 'MinPeakHeight', -abs(thresh), 'MinPeakDistance',200, 'MinPeakProminence',15);

obj2_app_DLC = camera_ts_FOR(obj2_app_DLC_frames);
obj1_app_DLC = camera_ts_FOR(obj1_app_DLC_frames);
figure
plot(camera_ts_FOR, -snout_obj2, obj2_app_DLC, dist_approach_obj2, 'o')
hold on 
%plot(obj2_app, -snout_obj2(events((events(:,3) == 2), 1) - events(2,1)), 'x','MarkerSize',15)
legend('distance from obj2', 'Auto DLC', 'Manual')
saveas(gcf, fullfile(path, 'plots_DLC', sprintf('Object_2_approach.png')))
figure
plot(camera_ts_FOR, -snout_obj1, obj1_app_DLC, dist_approach_obj1, 'o')
hold on
%plot(obj1_app, -snout_obj1(events(events(:,3) == 1, 1) - events(2,1)), 'x', 'MarkerSize',15)
legend('distance from obj1', 'Auto DLC', 'Manual')
saveas(gcf, fullfile(path, 'plots_DLC', sprintf('Object_1_approach.png')))
writematrix(obj1_app_DLC_frames+events(2,1), fullfile(path, sprintf('Object_1_approach.csv')))
writematrix(obj2_app_DLC_frames+events(2,1), fullfile(path, sprintf('Object_2_approach.csv')))
