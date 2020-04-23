close all;

global ax1 p dims target theta changed mode xPID yPID xAPID yAPID kp1Field ki1Field kd1Field kp2Field ki2Field kd2Field;

% PHYSICAL PARAMETERS/CONSTANTS
dims = [0.2, 0.2];
g = 9.8;
angleLim = 15;
maxSpeed = 180;

% INITIAL CONTROLLER SETUP
% THEORETICAL: kp1 = 462.6, kd1 = 257
% FAST: kp1 = 924, 129
% OPTIMIZED: kp1 = 763, kd1 = 257
% HAND-TUNED: kp1 = 800, kd1 = 500
kp1 = 763;
ki1 = 0;
kd1 = 257;
xPID = PIDController(kp1, ki1, kd1);
yPID = PIDController(kp1, ki1, kd1);
% THEORETICAL: kp2 = 28.6, kd2 = 0
% OPTIMIZED: kp2 = 28.6, kd2 = 0
% FAST: kp2 = 28.6, kd2 = 0
% HAND-TUNED: kp2 = 10, kd2 = 0
kp2 = 28.6;
ki2 = 0;
kd2 = 0;
xAPID = PIDController(kp2, ki2, kd2);
yAPID = PIDController(kp2, ki2, kd2);

% SIMULATION VIEW SETUP
fig1 = figure('Units', 'Normalized', 'Position', [0.1, 0.1, 0.8, 0.8], 'WindowButtonDownFcn', @onClick);
subplot(2, 2, 1);
hold on;
ax1 = gca;
ax1.Position(1) = -0.05;
ax1.Position(3) = 0.6;
ax1.Position(2) = 0.25;
ax1.Position(4) = 0.9;
axis off;
xlim([-0.3, 0.3]);
ylim([-0.3, 0.3]);
zlim([-0.3, 0.3]);
view([40, 20]);

% AXES SETUP FOR PLOTS
ax2 = subplot(2, 2, 2);
hold on;
grid on;
title('Position Controllers');
set(ax2, 'fontsize', 14);
ax3 = subplot(2, 2, 4);
hold on;
grid on;
title('Angle Controllers');
set(ax3, 'fontsize', 14);

% SETTINGS SETUP
bg = uibuttongroup(gcf, 'units', 'normalized', 'position', [0.025, 0.05, 0.15, 0.25]);
uicontrol(bg, 'style', 'text', 'String', 'Target Selection', 'units', 'normalized', 'fontsize', 14, 'position', [0.05, 0.8, 0.9, 0.2]);
b1 = uicontrol(bg, 'style', 'radiobutton', 'String', 'Manual', 'units', 'normalized', 'fontsize', 14, 'userdata', 0, 'position', [0.05, 0.65, 0.9, 0.2], 'callback', @radio);
b3 = uicontrol(bg, 'style', 'radiobutton', 'String', 'Square', 'units', 'normalized', 'fontsize', 14, 'userdata', 1, 'position', [0.05, 0.5, 0.9, 0.2], 'callback', @radio);
b4 = uicontrol(bg, 'style', 'radiobutton', 'String', 'Circle', 'units', 'normalized', 'fontsize', 14, 'userdata', 2, 'position', [0.05, 0.35, 0.9, 0.2], 'callback', @radio);
b5 = uicontrol(bg, 'style', 'radiobutton', 'String', 'Random', 'units', 'normalized', 'fontsize', 14, 'userdata', 3, 'position', [0.05, 0.2, 0.9, 0.2], 'callback', @radio);
b6 = uicontrol(bg, 'style', 'pushbutton', 'String', 'Reset', 'units', 'normalized', 'fontsize', 14, 'userdata', 4, 'position', [0.05, 0.025, 0.9, 0.2], 'callback', @radio);

panel1 = uipanel(gcf, 'units', 'normalized', 'position', [0.2, 0.05, 0.3, 0.25], 'fontsize', 14);
uicontrol(panel1, 'style', 'text', 'String', 'Controller Gains', 'units', 'normalized', 'fontsize', 14, 'position', [0.05, 0.8, 0.9, 0.2]);
panel2 = uipanel(panel1, 'units', 'normalized', 'position', [0.025, 0.025, 0.3, 0.75]);
uicontrol(panel2, 'style', 'text', 'String', 'Position', 'units', 'normalized', 'fontsize', 14, 'position', [0.05, 0.8, 0.9, 0.2]);
uicontrol(panel2, 'style', 'text', 'String', 'Kp', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.55, 0.3, 0.2]);
uicontrol(panel2, 'style', 'text', 'String', 'Ki', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.3, 0.3, 0.2]);
uicontrol(panel2, 'style', 'text', 'String', 'Kd', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.05, 0.3, 0.2]);
kp1Field = uicontrol(panel2, 'style', 'edit', 'string', kp1, 'fontsize', 14, 'userdata', 1, 'units', 'normalized', 'position', [0.3, 0.55, 0.65, 0.2], 'callback', @gains);
ki1Field = uicontrol(panel2, 'style', 'edit', 'string', ki1, 'fontsize', 14, 'userdata', 2, 'units', 'normalized', 'position', [0.3, 0.3, 0.65, 0.2], 'callback', @gains);
kd1Field = uicontrol(panel2, 'style', 'edit', 'string', kd1, 'fontsize', 14, 'userdata', 3, 'units', 'normalized', 'position', [0.3, 0.05, 0.65, 0.2], 'callback', @gains);

panel3 = uipanel(panel1, 'units', 'normalized', 'position', [0.35, 0.025, 0.3, 0.75]);
uicontrol(panel3, 'style', 'text', 'String', 'Angle', 'units', 'normalized', 'fontsize', 14, 'position', [0.05, 0.8, 0.9, 0.2]);
uicontrol(panel3, 'style', 'text', 'String', 'Kp', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.55, 0.3, 0.2]);
uicontrol(panel3, 'style', 'text', 'String', 'Ki', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.3, 0.3, 0.2]);
uicontrol(panel3, 'style', 'text', 'String', 'Kd', 'units', 'normalized', 'fontsize', 14, 'position', [0.025, 0.05, 0.3, 0.2]);
kp2Field = uicontrol(panel3, 'style', 'edit', 'string', kp2, 'fontsize', 14, 'userdata', 4, 'units', 'normalized', 'position', [0.3, 0.55, 0.65, 0.2], 'callback', @gains);
ki2Field = uicontrol(panel3, 'style', 'edit', 'string', ki2, 'fontsize', 14, 'userdata', 5, 'units', 'normalized', 'position', [0.3, 0.3, 0.65, 0.2], 'callback', @gains);
kd2Field = uicontrol(panel3, 'style', 'edit', 'string', kd2, 'fontsize', 14, 'userdata', 6, 'units', 'normalized', 'position', [0.3, 0.05, 0.65, 0.2], 'callback', @gains);

panel4 = uipanel(panel1, 'units', 'normalized', 'position', [0.675, 0.025, 0.3, 0.75]);
uicontrol(panel4, 'style', 'text', 'String', 'Presets', 'units', 'normalized', 'fontsize', 14, 'position', [0.05, 0.8, 0.9, 0.2]);
uicontrol(panel4, 'style', 'pushbutton', 'string', 'Optimized', 'fontsize', 14, 'userdata', 0, 'units', 'normalized', 'position', [0.05, 0.625, 0.9, 0.2], 'callback', @gains);
uicontrol(panel4, 'style', 'pushbutton', 'string', 'Theoretical', 'fontsize', 14, 'userdata', 0, 'units', 'normalized', 'position', [0.05, 0.425, 0.9, 0.2], 'callback', @gains);
uicontrol(panel4, 'style', 'pushbutton', 'string', 'Fast', 'fontsize', 14, 'userdata', 0, 'units', 'normalized', 'position', [0.05, 0.225, 0.9, 0.2], 'callback', @gains);
uicontrol(panel4, 'style', 'pushbutton', 'string', 'Hand-Tuned', 'fontsize', 14, 'userdata', 0, 'units', 'normalized', 'position', [0.05, 0.025, 0.9, 0.2], 'callback', @gains);

% POSITION PLOT SETUP
xPlot = animatedline(ax2, 0, 0, 'color', 'b', 'linestyle', '-', 'linewidth', 2);
xTargetPlot = animatedline(ax2, [0, 0], [0, 0], 'color', 'b', 'linestyle', ':', 'linewidth', 2);
yPlot = animatedline(ax2, 0, 0, 'color', 'r', 'linestyle', '-', 'linewidth', 2);
yTargetPlot = animatedline(ax2, [0, 0], [0, 0], 'color', 'r', 'linestyle', ':', 'linewidth', 2);
xlim(ax2, [0, 10]);
ylim(ax2, [-0.25, 0.25]);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Position (m)');

% ANGLE PLOT SETUP
xAnglePlot = animatedline(ax3, 0, 0, 'color', 'b', 'linestyle', '-', 'linewidth', 2);
xAngleTargetPlot = animatedline(ax3, [0, 0], [0, 0], 'color', 'b', 'linestyle', ':', 'linewidth', 2);
yAnglePlot = animatedline(ax3, 0, 0, 'color', 'r', 'linestyle', '-', 'linewidth', 2);
yAngleTargetPlot = animatedline(ax3, [0, 0], [0, 0], 'color', 'r', 'linestyle', ':', 'linewidth', 2);
xlim(ax3, [0, 10]);
ylim(ax3, [-20, 20]);
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Angle (°)');

% INITIAL VALUES
theta = [0, 0];
dTheta = [0, 0];
ddTheta = [0, 0];
thetaDes = [0, 0];
pos = [0, 0];
target = [0, 0];
v = [0, 0];
t = 0;
dt = 0.01;
changed = false;
[p, legs] = drawPlatform(theta);
b = drawBall(pos, theta);
T = drawTarget(target, theta);
mode = 0;
squareCount = 1;
squareCorners = [0.8, 0.8; 0.8, -0.8; -0.8, -0.8; -0.8, 0.8];
octCount = 1;
octCorners = 0.8 * [cos(0:pi/16:2*pi-0.01)',sin(0:pi/16:2*pi-0.01)'];

% MAIN LOOP
tic;
while ishandle(1)
    t = t + dt;
    nextTime = t + dt;
    
    % Reset controllers if target has changed
    if changed
        xPID.reset();
        yPID.reset();
        xAPID.reset();
        yAPID.reset();
        changed = false;
    end
    
    % select target based on mode
    if mode == 1
        if sum(abs(pos - target)) < 0.001 && sum(abs(v)) < 0.01
            target = squareCorners(squareCount, :) .* dims;
            squareCount = mod(squareCount, 4) + 1;
        end
    elseif mode == 2
        currentAngle = atan2d(pos(2), pos(1));
        nextAngle = currentAngle + 2;
        target = 0.1 * [cosd(nextAngle), sind(nextAngle)];
    elseif mode == 3
        if sum(abs(pos - target)) < 0.001 && sum(abs(v)) < 0.001
            target = (rand(1, 2) .* dims * 2 - dims) * 0.9;
        end
    elseif mode == 4
        target = [0, 0];
        mode = 0;
        b1.Value = 1;
    end
    
    % Get desired angles from position controllers
    e1 = pos - target;
    thetaDes(1) = xPID.control(t, e1(1));
    thetaDes(2) = yPID.control(t, e1(2));
    thetaDes(abs(thetaDes) > angleLim) = sign(thetaDes(abs(thetaDes) > angleLim)) * angleLim;
    
    % Get rotation speeds from angle controllers
    e2 = (thetaDes - theta) * pi / 180;
    dTheta(1) = xAPID.control(t, e2(1));
    dTheta(2) = yAPID.control(t, e2(2));
    dTheta = dTheta * 180 / pi;
    dTheta(dTheta > maxSpeed) = maxSpeed;
    dTheta(dTheta < -maxSpeed) = -maxSpeed;
    
    % Turn platform
    theta = theta + dTheta * dt;
    theta(abs(theta) > angleLim) = sign(theta(abs(theta) > angleLim)) * angleLim;
    
    % Move ball
    a = -5/7 * g * sind(theta);
    v = v + a * dt;    
    pos = pos + v * dt;
    
    % Redraw platform, target, and ball if anything has moved
    if any(v ~= 0) || any(dTheta ~= 0)
        [p, legs] = drawPlatform(theta, p, legs);
        T = drawTarget(target, theta, T);
        b = drawBall(pos, theta, b);
    end
    
    % Plot new ball position, target position, platform angles, and target angles
    addpoints(xPlot, t, pos(1));
    addpoints(yPlot, t, pos(2));
    addpoints(xTargetPlot, t, target(1));
    addpoints(yTargetPlot, t, target(2));
    addpoints(xAnglePlot, t, theta(1));
    addpoints(yAnglePlot, t, theta(2));
    addpoints(xAngleTargetPlot, t, thetaDes(1));
    addpoints(yAngleTargetPlot, t, thetaDes(2));
    
    % Pan plots
    if t > 9
        xlim(ax2, [t - 9, t + 1]);
        xlim(ax3, [t - 9, t + 1]);
    end
    
    % Pause
    waitUntil(nextTime)
end

% DRAW PLATFORM
function [p, legs] = drawPlatform(theta, p, legs)
global ax1 dims;

% Find corner locations
W = dims(1);
H = dims(2);
x = cosd(theta(1)) * [-W, W; -W, W];
y = cosd(theta(2)) * [-H, -H; H, H];
z = sind(theta(1)) * [-W, W; -W, W] + sind(theta(2)) * [-H, -H; H, H];

% Find leg endpoint positions
xl = [0.2, 0.2, -0.2; x(3), x(4), (x(1) + x(2)) / 2] * 0.9;
yl = [-0.2, 0.2, 0; y(3), y(4), (y(1) + y(2)) / 2] * 0.9;
zl = [-0.2, -0.2, -0.2; z(3), z(4), (z(1) + z(2)) / 2] * 0.9;
zl = (zl + 0.2) * 0.99 - 0.2;

% Get legs as vectors
leg1 = [xl(2) - xl(1), yl(2) - yl(1), zl(2)- zl(1)];
leg2 = [xl(4) - xl(3), yl(4) - yl(3), zl(4)- zl(3)];
leg3 = [xl(6) - xl(5), yl(6) - yl(5), zl(6)- zl(5)];

% Get perpendicular vector at the base of each leg.
base1 = [-1, -1, 0];
base2 = [1, -1, 0];
base3 = [0, 1, 0];

% Find the length of each leg
len1 = sqrt((xl(1) - xl(2)) ^ 2 + (yl(1) - yl(2)) ^ 2 + (zl(1) - zl(2)) ^ 2);
len2 = sqrt((xl(3) - xl(4)) ^ 2 + (yl(3) - yl(4)) ^ 2 + (zl(3) - zl(4)) ^ 2);
len3 = sqrt((xl(5) - xl(6)) ^ 2 + (yl(5) - yl(6)) ^ 2 + (zl(5) - zl(6)) ^ 2);

% Determine how far the middle joint will stick out
lengths = 0.14;
disp1 = lengths * sin(acos(len1 / (2 * lengths)));
disp2 = lengths * sin(acos(len2 / (2 * lengths)));
disp3 = lengths * sin(acos(len3 / (2 * lengths)));

% Find vector perpendicular to the leg and the base
orth1 = cross(leg1, base1);
orth2 = cross(leg2, base2);
orth3 = cross(leg3, base3);
orth1 = disp1 * orth1 / norm(orth1);
orth2 = disp2 * orth2 / norm(orth2);
orth3 = disp3 * orth3 / norm(orth3);

% Find location of the middle joint of each leg
midpt1 = [xl(1) + leg1(1) / 2 + orth1(1), yl(1) + leg1(2) / 2 + orth1(2), zl(1) + leg1(3) / 2 + orth1(3)];
midpt2 = [xl(3) + leg2(1) / 2 + orth2(1), yl(3) + leg2(2) / 2 + orth2(2), zl(3) + leg2(3) / 2 + orth2(3)];
midpt3 = [xl(5) + leg3(1) / 2 + orth3(1), yl(5) + leg3(2) / 2 + orth3(2), zl(5) + leg3(3) / 2 + orth3(3)];

% Find the x, y, and z values of all leg joints
xl = [xl(1, :); midpt1(1), midpt2(1), midpt3(1); xl(2, :)];
yl = [yl(1, :); midpt1(2), midpt2(2), midpt3(2); yl(2, :)];
zl = [zl(1, :); midpt1(3), midpt2(3), midpt3(3); zl(2, :)];

% if the platform and legs have already been created, move them. If not, create them.
if nargin == 3
    p.XData = x;
    p.YData = y;
    p.ZData = z;
    legs(1).XData = xl(:, 1);
    legs(2).XData = xl(:, 2);
    legs(3).XData = xl(:, 3);
    legs(1).YData = yl(:, 1);
    legs(2).YData = yl(:, 2);
    legs(3).YData = yl(:, 3);
    legs(1).ZData = zl(:, 1);
    legs(2).ZData = zl(:, 2);
    legs(3).ZData = zl(:, 3);
else
    p = surf(ax1, x, y, z, 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 1);
    legs = plot3(ax1, xl, yl, zl, 'k-', 'linewidth', 20);
end
end

% DRAW TARGET
function T = drawTarget(target, theta, T)
global ax1;

% Find the exact location of the target circle
target = target .* cosd(theta);
z = sind(theta(1)) * target(1) + sind(theta(2)) * target(2) + 0.018;

% If the target has already been created, move it. If not, create it.
if nargin == 3
    T.XData = target(1);
    T.YData = target(2);
    T.ZData = z;
else
    T = plot3(ax1, target(1), target(2), z, 'go', 'MarkerSize', 30, 'linewidth', 3);
end
end

% DRAW BALL
function b = drawBall(pos, theta, b)
global ax1;

% Find the exact location of the ball
pos3d = [pos .* cosd(theta), sum(pos .* sind(theta)) + 0.018];

% If the ball has already been created, move it. If not, create it.
if nargin == 3
    b.XData = pos3d(1);
    b.YData = pos3d(2);
    b.ZData = pos3d(3);
else
    b = plot3(ax1, pos3d(1), pos3d(2), pos3d(3), 'k.', 'MarkerSize', 100);
end
end

% CALLBACK FUNCTION FOR CLICKING ON THE PLATFORM
function onClick(~, ~)
global ax1 p dims target theta changed mode;

% Ensure that it is in manual target mode
if mode ~= 0
    return;
end

% Get vector representing click
line = get(ax1, 'CurrentPoint');
p1 = [p.XData(1), p.YData(1), p.ZData(1)];
p2 = [p.XData(2), p.YData(2), p.ZData(2)];
p3 = [p.XData(3), p.YData(3), p.ZData(3)];
p4 = line(1, :);
p5 = line(2, :);

% Find intersection between vector and plate
normal = cross(p1 - p2, p1 - p3);
syms x y z t
plane = dot(normal, [x, y, z] - p1);
line2 = p4 + t * (p5 - p4);
fcn = subs(plane, [x, y, z], line2);
t0 = solve(fcn);
point = subs(line2, t, t0);
point = double(point(1:2));
target2 = point ./ cosd(theta);

% Ensure the clicked point is on the plate
if abs(target2) < dims
    target = target2;
    changed = true;
end
end

% WAIT FUNCTION
function waitUntil(time)
% Wait until the time has been reached. (More accurate than simply using pause).
while toc < time
    pause(0.001);
end
end

% CALLBACK FOR MODE RADIO BUTTONS
function radio(~, event)
global mode;

mode = event.Source.UserData;
end

% CALLBACK FOR GAIN INPUT FIELDS
function gains(~, event)
global xPID yPID xAPID yAPID kp1Field ki1Field kd1Field kp2Field ki2Field kd2Field;

% Check if a preset button was pressed. If so, set all gains to their appropriate values.
if event.Source.UserData == 0
    switch event.Source.String
        case 'Optimized'
            kp1 = 763;
            ki1 = 0;
            kd1 = 257;
            kp2 = 28.6;
            ki2 = 0;
            kd2 = 0;
            kp1Field.String = kp1;
            ki1Field.String = ki1;
            kd1Field.String = kd1;
            kp2Field.String = kp2;
            ki2Field.String = ki2;
            kd2Field.String = kd2;
            xPID.setGains(kp1, ki1, kd1);
            yPID.setGains(kp1, ki1, kd1);
            xAPID.setGains(kp2, ki2, kd2);
            yAPID.setGains(kp2, ki2, kd2);
        case 'Theoretical'
            kp1 = 46;
            ki1 = 0;
            kd1 = 255;
            kp2 = 15.9;
            ki2 = 0;
            kd2 = 0;
            kp1Field.String = kp1;
            ki1Field.String = ki1;
            kd1Field.String = kd1;
            kp2Field.String = kp2;
            ki2Field.String = ki2;
            kd2Field.String = kd2;
            xPID.setGains(kp1, ki1, kd1);
            yPID.setGains(kp1, ki1, kd1);
            xAPID.setGains(kp2, ki2, kd2);
            yAPID.setGains(kp2, ki2, kd2);
        case 'Fast'
            kp1 = 924;
            ki1 = 0;
            kd1 = 129;
            kp2 = 157;
            ki2 = 0;
            kd2 = 0;
            kp1Field.String = kp1;
            ki1Field.String = ki1;
            kd1Field.String = kd1;
            kp2Field.String = kp2;
            ki2Field.String = ki2;
            kd2Field.String = kd2;
            xPID.setGains(kp1, ki1, kd1);
            yPID.setGains(kp1, ki1, kd1);
            xAPID.setGains(kp2, ki2, kd2);
            yAPID.setGains(kp2, ki2, kd2);
        case 'Hand-Tuned'
            kp1 = 800;
            ki1 = 0;
            kd1 = 500;
            kp2 = 10;
            ki2 = 0;
            kd2 = 0.5;
            kp1Field.String = kp1;
            ki1Field.String = ki1;
            kd1Field.String = kd1;
            kp2Field.String = kp2;
            ki2Field.String = ki2;
            kd2Field.String = kd2;
            xPID.setGains(kp1, ki1, kd1);
            yPID.setGains(kp1, ki1, kd1);
            xAPID.setGains(kp2, ki2, kd2);
            yAPID.setGains(kp2, ki2, kd2);
    end
    return;
end

% Get number inputted in input field, and exit if it was not a number.
num = str2double(event.Source.String);
if isempty(num)
    event.Source.String = 0;
    return;
end

% Set the appropriate gain in the appropriate controller.
switch event.Source.UserData
    case 1
        xPID.setGains('p', num);
        yPID.setGains('p', num);
    case 2
        xPID.setGains('i', num);
        yPID.setGains('i', num);
    case 3
        xPID.setGains('d', num);
        yPID.setGains('d', num);
    case 4
        xAPID.setGains('p', num);
        yAPID.setGains('p', num);
    case 5
        xAPID.setGains('i', num);
        yAPID.setGains('i', num);
    case 6
        xAPID.setGains('d', num);
        yAPID.setGains('d', num);
end
end