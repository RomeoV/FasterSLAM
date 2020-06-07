clear; clc; close all;

% Number of Waypoints ( Corners )
nwp = 25;
% Noise level one x and y coordinates ( controls the distortion of the reference circle )
radius = 100;
x_noise = 70*rand(nwp,1);
y_noise = 70*rand(nwp,1);

wp = create_wp(0, 0, radius, nwp, x_noise, y_noise);
lm = create_lm(-(radius+radius/2), (radius+radius/2), 200, wp);
% Specify seed
%seed = 10;
%rng(seed);

figure, hold on; grid;
plot([wp(:,1);wp(1,1)], [wp(:,2);wp(1,2)], 'k*');
plot([wp(:,1);wp(1,1)], [wp(:,2);wp(1,2)], 'r');
plot(lm(:,1), lm(:,2), 'mo');
hold off;
write_file(lm, wp, 'input_4.txt');

function wp = create_wp(x_center, y_center, r, nwp, x_noise, y_noise)
    % Waypoints
    step_size = 2*pi / nwp;
    theta = ( 0 : step_size : 2*pi - step_size )';
    x = r * cos(theta) + x_center + x_noise;
    y = r * sin(theta) + y_center + y_noise; 
    wp = [x y];
    wp = wp - mean(wp);
end

function lm = create_lm(a, b, nlm, wp)
    lm = rand(nlm,2)*(b-a) + a;
    mask = false(length(lm),1);
    for i = 1:length(lm)
        dist = sqrt( sum( ( wp - lm(i,:) ).^2, 2 ) );
        % Landmark thresholding
        if min(dist) > 25 && min(dist) < 50
            mask(i) = true;
        end
    end 
    lm = lm(mask, :);
end

function write_file(lm, wp, fname)
    flag = input(['Should I write the data in ' fname ' ? [1/0] ']);
    if flag ~= 1
        return;
    end
    fid = fopen(fname, 'w');
    fprintf(fid, '#type rows columns\n');
    fprintf(fid, 'lm %d %d\n', size(lm,1), size(lm,2));
    for i = 1 : size(lm,1)
        fprintf(fid, '%f %f\n', lm(i,1), lm(i,2));
    end
    fprintf(fid, '\nwp %d %d\n', size(wp,1), size(wp,2));
    for i = 1 : size(wp,1)
        fprintf(fid, '%f %f\n', wp(i,1), wp(i,2));
    end
    fclose(fid);
end
