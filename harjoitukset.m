% Algebra ja geometria (4op)
% Octave-harjoitukset

% Trigonometria: 1-8, 11, 15
% Analyyttinen geometria: 1-3, 7, 8
% Analyyttinen geometria osa 2: 1-4, 6-9
% Sinikäyrät: 1-3
% Kompleksiluvut: 1
% eksponentti ja logaritmi: 1-2

% yhteensä 29



%trig_1

clear
L = 6
x = 2
dx = 1

y1 = sqrt(L^2-x^2);
y2 = sqrt(L^2 - (x + dx)^2);

dy = y1-y2



%trig_2

clear
L = 5;
x = 7;

OP = sqrt(L^2 - (x/2)^2)
OP = OP * 2



%trig_3

clear
h = 1.5;
L = 6;

R = (h^2 + (L^2)/4) / (2*h);

disp(R);



%trig_4

clear
R = 5;   % Radius of the sphere
r = 2.5; % Radius of the cylinder base
h = 20;  % Height of the cylinder

% Calculate the height of the spherical cap (k)
k = sqrt(R^2 - r^2);

% Cylinder surface area
lpA = pi * r^2; % Base area of the cylinder
lk = h - 2 * k; % Height of the cylindrical portion excluding the sphere
lsA = 2 * pi * r * lk; % Lateral surface area of the cylinder
lA = lsA + 2 * lpA; % Total surface area of the cylinder including both bases

% Sphere surface area
palloA = 4 * pi * R^2; % Full surface area of the sphere
palloSA = 2 * pi * R * (R - k); % Surface area of two spherical caps
A = (palloA - 2 * palloSA) + lA; % Total surface area

% Volume calculations
pV = (4 / 3) * pi * R^3; % Volume of the sphere
lV1 = pi * r^2 * lk; % Volume of the cylindrical part

% Volume of the spherical caps
h_cap = R - k; % Height of one spherical cap
segmenttiV = 2 * (pi / 3) * h_cap^2 * (3 * R - h_cap); % Volume of two spherical caps
V = pV + lV1 - segmenttiV; % Total volume



%trig_5

clear
h = 2;
H = 3;
L = 8;
v1 = 1;
v2 = 2;
x = 0:0.01:L;

AP = sqrt(x.^2 + h^2);
PB = sqrt(H^2 + (L - x).^2);

APt = AP / v1;
PBt = PB / v2;
tt = APt + PBt;

% Find the minimum time and corresponding x value
[min_time, min_index] = min(tt);
min_x = x(min_index);

% Configures the plot
figure(1)
plot(x, tt, 'b', 'LineWidth', 1)
hold on
plot(min_x, min_time, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r')
ylabel('aika t')
xlabel('x')
grid on
title(sprintf('h = %d, H = %d, L = %d, v_1 = %d, v_2: = %d tmin = 6.04, kun x = 1.03', h, H, L, v1, v2))

% Set ticks
xticks(0:1:L)
yticks(6:0.5:10)


hold off



%trig_6

clear
r = 1
R = 3
a = 220

r1 = (a/360)*r
r2 = (a/360)*R

h = sqrt((R-r)^2 - (r2-r1)^2)



%trig_7

clear
% Given values
h = 4;
alpha = 35;
beta = 10;

% Convert degrees to radians
alpha_rad = deg2rad(alpha);
beta_rad = deg2rad(beta);

% Calculate x
x = h * tan(beta_rad);

% Calculate L
L = (h * tan(alpha_rad+beta_rad)) - x;

% Display results
fprintf('x = %.3f\n', x);
fprintf('L = %.3f\n', L);



%trig_8

clear
% Given values
a = 1;
b = 2;
p = 0.2;
alpha_values = linspace(10, 80, 100);  % angles from 10 to 80 degrees
alpha_rad = deg2rad(alpha_values);     % convert degrees to radians

% Calculate AB, AP, and QB based on the given formula
AB = b ./ sin(alpha_rad) + a ./ cos(alpha_rad);
AP = p ./ sin(alpha_rad);
QB = p ./ cos(alpha_rad);

% CD is computed as
CD = AB - AP - QB;

% Find the minimum value of CD and the corresponding alpha
[CDmin, idx_min] = min(CD);
alpha_min = alpha_values(idx_min);

% Plotting
figure;
plot(alpha_values, CD, 'b', 'LineWidth', 2);
hold on;
plot(alpha_min, CDmin, 'ro', 'MarkerSize', 8);  % Mark the minimum point

title('a = 1, b = 2, p = 0.2, CDmin = 3.7495, kun \alpha = 52.6^\circ');
xlabel('kulma \alpha');
ylabel('CD');
grid on;

% Display the minimum value of CD and corresponding alpha
disp(['Minimum CD: ', num2str(CDmin), ' at alpha = ', num2str(alpha_min), ' degrees']);



%trig_11

clear
a = 7;
b = 4;
c = 5;
d = 3;
h = 4;
z = 1;

e = c + ((a - c) / h) * z;
f = d + ((b - d) / h) * z;

V_neste = (1/6) * (2*e*f + e*d + f*c + 2*c*d) * z;

V_neste



%trig_15

clear
function common_area = overlap_area(R, r, L)
  if L >= R + r
    common_area = 0; % No overlap
  elseif L <= abs(R - r)
    common_area = pi * min(R, r)^2; % One circle is completely inside the other
  else
    % Calculate the area of overlap
    part1 = R^2 * acos((L^2 + R^2 - r^2) / (2 * L * R));
    part2 = r^2 * acos((L^2 + r^2 - R^2) / (2 * L * r));
    part3 = 0.5 * sqrt((-L + R + r) * (L + R - r) * (L - R + r) * (L + R + r));
    common_area = part1 + part2 - part3;
  end
end

% Example usage
R = 5;
r = 3;
L = 6;
common_area = overlap_area(R, r, L);
disp(common_area);



%anageo_1

clear
A = [1, 1];
B = [5, 2];
C = [4, 4];

% Lengths of the sides
AB = norm(B - A);
BC = norm(C - B);
CA = norm(A - C);

% Semi-perimeter of the triangle
p = (AB + BC + CA) / 2;

% Area of the triangle using Heron's formula
area = sqrt(p * (p - AB) * (p - BC) * (p - CA));

% Calculating angles using the cosine rule
alpha = acosd((BC^2 + CA^2 - AB^2) / (2 * BC * CA));
beta = acosd((AB^2 + CA^2 - BC^2) / (2 * AB * CA));
gamma = acosd((AB^2 + BC^2 - CA^2) / (2 * AB * BC));

% Display the results
fprintf('AB = %.4f, BC = %.4f, CA = %.4f\n', AB, BC, CA);
fprintf('alpha = %.4f°, beta = %.4f°, gamma = %.4f°\n', alpha, beta, gamma);
fprintf('Area = %.4f\n', area);

% Plot the triangle and points

% Plot the triangle edges
plot([A(1), B(1)], [A(2), B(2)], 'b-', 'LineWidth', 2); % AB
hold on;
plot([B(1), C(1)], [B(2), C(2)], 'b-', 'LineWidth', 2); % BC
plot([C(1), A(1)], [C(2), A(2)], 'b-', 'LineWidth', 2); % CA

% Plot the points A, B, C as scatter points
hA = scatter(A(1), A(2), 50, 'r', 'filled'); % Point A
hB = scatter(B(1), B(2), 50, 'g', 'filled'); % Point B
hC = scatter(C(1), C(2), 50, 'b', 'filled'); % Point C

% Create a legend specifically for the points
legend([hA, hB, hC], {'A', 'B', 'C'}, 'Location', 'best');

% Add grid and formatting
grid on;
xticks(0:0.5:5);
yticks(0:0.5:5);
title(sprintf('AB = %.4f, BC = %.4f, CA = %.4f\n\\alpha = %.4f°, \\beta = %.4f°, \\gamma = %.4f°, Area = %.4f', AB, BC, CA, beta, gamma, alpha, area));



%anageo_2

clear
% Coordinates of points A, B, and C
A = [1, 1];
B = [5, 2];
C = [3, 5];

% Function to calculate the slope between two points
slope = @(p1, p2) (p2(2) - p1(2)) / (p2(1) - p1(1));

% Slopes of sides BC, CA, and AB
m_BC = slope(B, C);
m_CA = slope(C, A);
m_AB = slope(A, B);

% Slopes of altitudes (perpendicular slopes)
m_AH = -1 / m_BC;
m_BI = -1 / m_CA;
m_CJ = -1 / m_AB;

% Calculate intercepts of altitudes using point-slope form y - y1 = m(x - x1)
b_AH = A(2) - m_AH * A(1);
b_BI = B(2) - m_BI * B(1);
b_CJ = C(2) - m_CJ * C(1);

% Finding the intersection points H, I, J (feet of altitudes)
% Intersection of altitude AH with BC
H_x = (b_AH - (C(2) - m_BC * C(1))) / (m_BC - m_AH);
H_y = m_AH * H_x + b_AH;
H = [H_x, H_y];

% Intersection of altitude BI with CA
I_x = (b_BI - (A(2) - m_CA * A(1))) / (m_CA - m_BI);
I_y = m_BI * I_x + b_BI;
I = [I_x, I_y];

% Intersection of altitude CJ with AB
J_x = (b_CJ - (A(2) - m_AB * A(1))) / (m_AB - m_CJ);
J_y = m_CJ * J_x + b_CJ;
J = [J_x, J_y];

% Finding orthocenter P (intersection of two altitudes AH and BI)
P_x = (b_AH - b_BI) / (m_BI - m_AH);
P_y = m_AH * P_x + b_AH;
P = [P_x, P_y];

% Display the coordinates of H, I, J, and P
H, I, J, P

% Plotting the triangle ABC, altitudes, and orthocenter
figure;
hold on;
plot([A(1), B(1)], [A(2), B(2)], 'b-', 'LineWidth', 1);
plot([B(1), C(1)], [B(2), C(2)], 'b-', 'LineWidth', 1);
plot([C(1), A(1)], [C(2), A(2)], 'b-', 'LineWidth', 1);

% Plot altitudes as solid lines in red
plot([A(1), H(1)], [A(2), H(2)], 'r-', 'LineWidth', 1);
plot([B(1), I(1)], [B(2), I(2)], 'r-', 'LineWidth', 1);
plot([C(1), J(1)], [C(2), J(2)], 'r-', 'LineWidth', 1);

% Plot orthocenter P
plot(P(1), P(2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

grid on;
xticks(0:1:6);
yticks(0:1:6);
axis equal;
hold off;



%anageo_3

clear
% Define vertices of the triangle
A = [0, 0];
B = [5, 1];
C = [3, 6];

% Define sets of coordinates for P and Q
coordinate_pairs = {
    [1, 4], [5, 3];
    [2, 3.5], [5, 3];
    [2, 3.5], [4, 3];
    [1, 4], [-1, 5]
};

% Function to find intersection point
function [t, S, isValid] = line_intersect(P, Q, X, Y)
    dPQ = Q - P;
    dXY = Y - X;

    denom = dPQ(1) * dXY(2) - dPQ(2) * dXY(1);
    isValid = true;

    if abs(denom) < 1e-10
        t = NaN;  % Parallel or collinear
        S = NaN(1, 2);
        isValid = false; % No valid intersection if parallel
    else
        t = ((X(1) - P(1)) * dXY(2) - (X(2) - P(2)) * dXY(1)) / denom;
        S = P + t * dPQ;

        % Check if the intersection point S lies on the line segment X-Y
        % calculate parameter 'r' for segment XY
        r = ((S(1) - X(1)) * (Y(1) - X(1)) + (S(2) - X(2)) * (Y(2) - X(2))) / sum(dXY.^2);
        if r < 0 || r > 1
            isValid = false; % Intersection point is outside segment XY
        end
    end
end

% Function to check if a point is inside a triangle using barycentric coordinates
function inside = is_point_in_triangle(P, A, B, C)
    v0 = C - A;
    v1 = B - A;
    v2 = P - A;

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    inside = (u >= 0) && (v >= 0) && (u + v <= 1);
end

% Loop through each coordinate pair for P and Q
for i = 1:size(coordinate_pairs, 1)
    % Get P and Q for this iteration
    P = coordinate_pairs{i, 1};
    Q = coordinate_pairs{i, 2};

    % Check if either P or Q is inside the triangle
    P_inside = is_point_in_triangle(P, A, B, C);
    Q_inside = is_point_in_triangle(Q, A, B, C);

    % Check intersection with each triangle side
    [t1, S1, valid1] = line_intersect(P, Q, A, B);
    [t2, S2, valid2] = line_intersect(P, Q, B, C);
    [t3, S3, valid3] = line_intersect(P, Q, C, A);

    % Collect valid intersections within [0, 1] range for segment PQ
    intersections = [];
    if valid1 && t1 >= 0 && t1 <= 1, intersections = [intersections; S1]; end
    if valid2 && t2 >= 0 && t2 <= 1, intersections = [intersections; S2]; end
    if valid3 && t3 >= 0 && t3 <= 1, intersections = [intersections; S3]; end

    % Plotting
    figure;
    hold on;
    grid on
    xticks(-1:1:6);
    yticks(0:1:6);
    xlim([-1 6]);
    ylim([0 6]);
    h1 = plot([A(1) B(1) C(1) A(1)], [A(2) B(2) C(2) A(2)], 'b-', 'LineWidth', 1.5); % Triangle ABC

    % Determine how to color the PQ segment
    if P_inside && Q_inside
        % If both points are inside, plot the whole segment in normal green
        h2 = plot([P(1) Q(1)], [P(2) Q(2)], 'g-', 'LineWidth', 2); % Standard green
    elseif P_inside && ~isempty(intersections)
        % If only P is inside, draw from P to the first intersection in normal green
        h2 = plot([P(1) intersections(1, 1)], [P(2) intersections(1, 2)], 'g-', 'LineWidth', 2);
        plot([intersections(1, 1) Q(1)], [intersections(1, 2) Q(2)], 'r-', 'LineWidth', 1.5); % Remaining in red
    elseif Q_inside && ~isempty(intersections)
        % If only Q is inside, draw from Q to the first intersection in normal green
        h2 = plot([Q(1) intersections(1, 1)], [Q(2) intersections(1, 2)], 'g-', 'LineWidth', 2);
        plot([P(1) intersections(1, 1)], [P(2) intersections(1, 2)], 'r-', 'LineWidth', 1.5); % Remaining in red
    elseif ~P_inside && ~Q_inside && size(intersections, 1) == 2
        % Both endpoints are outside but there are two intersections
        % Plot the segments outside the triangle in red first
        plot([P(1) intersections(1, 1)], [P(2) intersections(1, 2)], 'r-', 'LineWidth', 1.5); % Before first intersection
        plot([Q(1) intersections(2, 1)], [Q(2) intersections(2, 2)], 'r-', 'LineWidth', 1.5); % After second intersection
        % Plot inside segment last in normal green to ensure visibility
        h2 = plot(intersections(:,1), intersections(:,2), 'g-', 'LineWidth', 2); % Inside in normal green
    else
        % Default to red if no part of PQ is inside the triangle
        h2 = plot([P(1) Q(1)], [P(2) Q(2)], 'r-', 'LineWidth', 1.5);
    end

    % Plot smaller circles for P and Q points
    h3 = plot(P(1), P(2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c', 'DisplayName', 'P'); % Smaller marker size
    h4 = plot(Q(1), Q(2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm', 'DisplayName', 'Q'); % Smaller marker size

    legend([h3, h4, h1], 'P', 'Q', 'ABC', 'Location', 'Best');
    hold off;
end



%anageo_7

clear;
% Coordinates of points A, B, and C
A = [1, 1];
B = [4, 4];
C = [5, 2];

% Function to calculate the distance between two points
distance = @(P, Q) norm(P - Q);

% Calculate side lengths AB, BC, and AC
AB = distance(A, B);
BC = distance(B, C);
AC = distance(A, C);

% Calculate x, y, z using the formulas
x = 0.5 * (AB + AC - BC);
y = 0.5 * (AB + BC - AC);
z = 0.5 * (BC + AC - AB);

% Calculate the incenter (I) as the weighted average of A, B, C
I = (A * BC + B * AC + C * AB) / (AB + BC + AC);

% Calculate the semiperimeter and the area of the triangle using Heron's formula
s = (AB + BC + AC) / 2;
area = sqrt(s * (s - AB) * (s - BC) * (s - AC));

% Calculate the radius of the incircle
r = area / s;

% Calculate points K, L, M (tangency points on each side of the triangle)
K = A + x * (B - A) / AB;  % point K on side AB
L = B + y * (C - B) / BC;  % point L on side BC
M = C + z * (A - C) / AC;  % point M on side AC

% Plot the triangle
figure;
hold on;
plot([A(1), B(1), C(1), A(1)], [A(2), B(2), C(2), A(2)], 'b-', 'LineWidth', 0.5);

% Plot the incircle
theta = linspace(0, 2*pi, 100);
incircle_x = I(1) + r * cos(theta);
incircle_y = I(2) + r * sin(theta);
plot(incircle_x, incircle_y, 'r-', 'LineWidth', 0.5);

% Plot the incenter and tangency points (K, L, M) only
plot(I(1), I(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);  % Incenter filled in red
plot([K(1), L(1), M(1)], [K(2), L(2), M(2)], 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3);  % Tangency points filled in blue

% Set plot limits and aspect ratio
axis equal;
xlim([0.5 5.5]);
ylim([0.5 4.5]);
xticks(0.5:0.5:5.5);
yticks(0.5:0.5:4.5);
grid on;

hold off;



%anageo_8

clear
% Define centers and radii of the circles
A = [1, 1]; % Center of first circle
B = [5, 3]; % Center of second circle
R = 4;      % Radius of first circle (red)
r = 3;      % Radius of second circle (blue)

% Calculate distance between centers
d = norm(B - A);

% Check if the circles intersect
if d > R + r
    intersect = false;
elseif d < abs(R - r)
    intersect = false;
elseif d == 0 && R == r
    intersect = false;
else
    % Circles intersect, calculate intersection points S and T
    intersect = true;
    % Distance from A to the midpoint P along line AB
    a = (R^2 - r^2 + d^2) / (2 * d);

    % Point P on line AB which is the midpoint of the intersection chord
    P = A + a * (B - A) / d;

    % Distance from P to the intersection points S and T
    h = sqrt(R^2 - a^2);

    % Calculate intersection points S and T
    offset = h * [-(B(2) - A(2)) / d, (B(1) - A(1)) / d];
    S = P + offset;
    T = P - offset;
end

% Plot the circles
theta = linspace(0, 2 * pi, 100);
circle_A_x = A(1) + R * cos(theta);
circle_A_y = A(2) + R * sin(theta);
circle_B_x = B(1) + r * cos(theta);
circle_B_y = B(2) + r * sin(theta);

figure;
hold on;
plot(circle_A_x, circle_A_y, 'r-', 'LineWidth', 1.5); % Circle A in red
plot(circle_B_x, circle_B_y, 'b-', 'LineWidth', 1.5); % Circle B in blue
plot(A(1), A(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3); % Center A
plot(B(1), B(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3); % Center B

if intersect
    % Plot the intersection points S and T with DisplayName for legend
    hS = plot(S(1), S(2), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 4); % S
    hT = plot(T(1), T(2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 4); % T

    % Add legend for S and T
    legend([hS, hT], 'S', 'T');
    title(['S = [', num2str(S(1), '%.2f'), ', ', num2str(S(2), '%.2f'), ...
           '], T = [', num2str(T(1), '%.2f'), ', ', num2str(T(2), '%.5f'), ']']);
else
    title('eivät leikkaa');
end

% Set plot limits and aspect ratio
axis equal;
xlim([-4, 9]);
ylim([-5, 8]);

% Set x and y ticks
xticks(-4:1:9);
yticks(-5:1:8);

grid on;
hold off;



%anageo2_1

clear
% Define parameters for the two lines
k = 3; b = 1;  % Line 1: y = kx + b
m = -2; c = 2; % Line 2: y = mx + c

% Solve for intersection point (x0, y0)
x0 = (c - b) / (k - m);
y0 = k * x0 + b;

% Generate x values for plotting
x = linspace(-2, 2, 100);

% Define the two lines
y1 = k * x + b;
y2 = m * x + c;

% Plot the two lines
figure;
plot(x, y1, 'r', 'LineWidth', 1.5); hold on; % Line 1 in red
plot(x, y2, 'b', 'LineWidth', 1.5);          % Line 2 in blue

% Highlight the intersection point
plot(x0, y0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % Intersection in black

% Add labels and legend
xlabel('x'); ylabel('y');
legend('y = kx + b', 'y = mx + c', '[x_0, y_0]', 'Location', 'SouthEast'); % Legend in bottom-right corner
title(sprintf('k = %d, b = %d, m = %d, c = %d: x_0 = %.1f, y_0 = %.1f', k, b, m, c, x0, y0));

% Set axis limits and grid
xticks(-2:0.5:5);
yticks(-6:2:8);
grid on;

hold off;



%anageo2_2

clear
x1 = 2
x2 = 6
y1 = 1
y2 = 7
z11 = 3.6
z12 = 1.9
z21 = 2.3
z22 = 1.4
x0 = 2.9
y0 = 3.2

function res = linear_interp1d(x, x0, x1, y0, y1)
  #syms x y x0 y0 x1 y1
  #solve((y - y0)/(x - x0) == (y1 - y0)/(x1 - x0), y)
  res = (x*y0 - x*y1 + x0*y1 - x1*y0)/(x0-x1);
end

z1 = linear_interp1d(x0, x1, x2, z11, z12)
z2 = linear_interp1d(x0, x1, x2, z21, z22)
z0 = linear_interp1d(y0, y1, y2, z1, z2)

subplot(1,3,1)
plot([x1, x1], [0, z11], "b-")
hold on
plot([x1, x2], [0, 0], "b-")
plot([x2, x2], [0, z12], "g-")
plot([x1, x2], [z11, z12], "k-")
plot([x0, x0], [0, z1], "r-")
axis([1 x2+1 0 z11])
axis square
grid
hold off

subplot(1,3,2)
plot([x1, x1], [0, z21], "m-")
hold on
plot([x1, x2], [0, 0], "b-")
plot([x2, x2], [0, z22], "c-")
plot([x1, x2], [z21, z22], "k-")
plot([x0, x0], [0, z2], "g-")
axis([1 x2+1 0 z21])
axis square
grid
hold off

subplot(1,3,3)
plot([y1, y1], [0, z1], "r-")
hold on
plot([y1, y2], [0, 0], "b-")
plot([y2, y2], [0, z2], "g-")
plot([y1, y2], [z1, z2], "b-")
plot([y0, y0], [0, z0], "k-")
axis([0 y2+1 0 z1])
axis square
grid
hold off

axes( 'visible', 'off', 'title', ["x_1 = " num2str(x1) ", x_2 = " num2str(x2) ", y_1 = "
num2str(y1) ", y_2 = " num2str(y2) "\nz_1_1 = " num2str(z22) ", z_1_2" num2str(z12)
", z_2_1 = " num2str(z21) ", z_2_2 = " num2str(z22) "\nx_0 = " num2str(x0) ", y_0 = "
num2str(y0) "\n-> z_1 = " num2str(z1) ", z_2 = " num2str(z2) ", z_0 = " num2str(z0)] );



%anageo2_3


clear
% Parameters
Px = 3;
Py = -2;

a = 0.2;
b = -0.5;
c = 0.5;

f1 = @(x) a * x.^2 + b*x + c; % Parabola
f2 = @(x, k) k * (x - Px) + Py; % Red lines

% Calculate slopes of tangent lines
k_values = [2*Px*a + b - 2*sqrt(a*(Px^2*a + Px*b - Py + c)), 2*Px*a + b + 2*sqrt(a*(Px^2*a + Px*b - Py + c))];
k_plus = k_values(1);
k_neg = k_values(2);

fprintf('k1 = %.3f, k2 = %.3f\n', k_values(1), k_values(2));

% Solve intersections
syms x
intersection1 = solve(f1(x) == f2(x, k_plus), x);
intersection2 = solve(f1(x) == f2(x, k_neg), x);

% Convert symbolic results to numeric
x_int1 = double(intersection1);
x_int2 = double(intersection2);

% Compute corresponding y values
y_int1 = double(f1(x_int1));
y_int2 = double(f1(x_int2));

% Plotting
figure;
plot(Px, Py, "r.", "MarkerSize", 15);
hold on
x_vals = linspace(-3, 9);
f1_plot = plot(x_vals, f1(x_vals), 'b', 'LineWidth', 1);
f2_plot = plot(x_vals, f2(x_vals, k_plus), 'r', 'LineWidth', 1);
f3_plot = plot(x_vals, f2(x_vals, k_neg), 'r', 'LineWidth', 1);

% Add black dots for intersections
plot(x_int1, y_int1, 'k.', 'MarkerSize', 10); % Black dot for intersection 1
plot(x_int2, y_int2, 'k.', 'MarkerSize', 10); % Black dot for intersection 2

hold off
grid
xticks(-5:1:11)
yticks(-2:1:11)
axis([-5 11 -2 11])
xlabel("x");
ylabel("y");
set(get(gca,'ylabel'),'rotation',0)

title(sprintf('P_x = %d, P_y = %d, a = %.1f, b = %.1f, c = %.1f\nk_1 = %.3f, k_2 = %.3f', ...
    Px, Py, a, b, c, k_values(1), k_values(2)));



%anageo2_4


clear
% Given data
t1 = 2; s1 = 4;  % First time and distance
t2 = 10; s2 = 12; % Second time and distance

% Rewrite equations into matrix form
% s1 = v0*t1 + (1/2)*a*t1^2
% s2 = v0*t2 + (1/2)*a*t2^2
% This becomes:
% [t1, 0.5*t1^2; t2, 0.5*t2^2] * [v0; a] = [s1; s2]

A = [t1, 0.5*t1^2; t2, 0.5*t2^2];
b = [s1; s2];

% Solve for [v0; a]
solution = A \ b;

% Extract v0 and a
v0 = solution(1);
a = solution(2);

% Display the results
fprintf('Correct values: v0 = %.2f, a = %.2f\n', v0, a);

% Define the distance function
t = linspace(0, t2+2, 1000); % Extend time range
s = v0*t + (1/2)*a*t.^2;

% Plot the distance vs time
figure;
plot(t, s, '-b', 'LineWidth', 2); hold on;
plot(t1, s1, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % Mark (t1, s1)
plot(t2, s2, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % Mark (t2, s2)

% Add labels, title, and legend
xlabel('Time (t)');
ylabel('Distance (s)');
title(sprintf('v_0 = %.2f, a = %.2f', v0, a));
grid on;

% Set axis limits
xlim([0, t2+2]);
ylim([0, max(s)+2]);



%anageo2_6

clear
% Parameters
L = 10; % Target horizontal coordinate
H = 6; % Target vertical coordinate
alpha_deg = 60; % Launch angle in degrees
g = 9.81; % Acceleration due to gravity (m/s^2)

% Convert angle to radians
alpha = deg2rad(alpha_deg);

% Correct calculation for minimum launch angle
min_alpha = atan(H / L); % Minimum angle in radians
fprintf("Minimum angle (min alpha): %.4f degrees\n", rad2deg(min_alpha));

% Check the minimum angle condition
if alpha < min_alpha
    error("The launch angle must be greater than %.4f degrees.", rad2deg(min_alpha));
end

% Calculate the initial velocity v0
v0 = sqrt((g * L^2) / (2 * (L * tan(alpha) - H) * cos(alpha)^2));

% Display the result
fprintf("Required initial velocity v0: %.4f m/s\n", v0);

% Coefficients for the trajectory equation
a = -g / (2 * v0^2 * cos(alpha)^2);
b = tan(alpha);

% Define the trajectory equation
x = linspace(0, L, 100); % Generate x values from 0 to L
y = a * x.^2 + b * x; % Compute corresponding y values

% Plot the trajectory
figure;
plot(x, y, 'b-', 'LineWidth', 2); % Trajectory
hold on;

% Add the target point
plot(L, H, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Target point

% Add the black line from the target to the ground
plot([L, L], [0, H], 'k-', 'LineWidth', 1); % Vertical line

% Add the black ground line
plot([0, L + 2], [0, 0], 'k-', 'LineWidth', 1); % Ground line

% Add labels, title, and legend
xlabel('x (m)');
ylabel('y (m)');
% Add a detailed title with mathematical expressions
% Add a detailed title with LaTeX formatting
title(['L = 10 , H = 6, a = 60, min alpha = 30.9638 : v_0 = 13.1649'], ...
      'Interpreter', 'latex', 'FontSize', 12);

grid on;
xlim([0, L + 2]);
ylim([0, max(y) + 2]);

hold off;



%anageo2_7

clear
% Define the three points
x1 = 1; y1 = 1;
x2 = 4; y2 = 4;
x3 = 5; y3 = 2;

% Solve for x0, y0 using two equations
A = [2*(x2-x1), 2*(y2-y1);
     2*(x3-x1), 2*(y3-y1)];
b = [(x2^2 + y2^2 - x1^2 - y1^2);
     (x3^2 + y3^2 - x1^2 - y1^2)];

center = A\b; % Solves the linear system
x0 = center(1);
y0 = center(2);

% Compute the radius
r = sqrt((x1 - x0)^2 + (y1 - y0)^2);

% Display the results
fprintf('Center: (%.2f, %.2f)\n', x0, y0);
fprintf('Radius: %.2f\n', r);

% Plot the circle and points
theta = linspace(0, 2*pi, 100);
x_circle = x0 + r * cos(theta);
y_circle = y0 + r * sin(theta);

figure;
plot(x_circle, y_circle, 'b'); % Circle
hold on;
plot([x1, x2, x3], [y1, y2, y3], 'ro'); % Points
plot(x0, y0, 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k'); % Center
axis equal;
grid on;
xlabel('x');
ylabel('y');



%anageo2_8

clear
% Given values
F = 4;          % Distance of focal points from the origin
Px = 2; Py = 3; % Point on the ellipse

% Focal points
F1 = [-F, 0];
F2 = [F, 0];

% Calculate distances
d1 = sqrt((Px - F1(1))^2 + (Py - F1(2))^2); % Distance to F1
d2 = sqrt((Px - F2(1))^2 + (Py - F2(2))^2); % Distance to F2

% Semi-major axis
a = (d1 + d2) / 2;

% Semi-minor axis
b = sqrt(a^2 - F^2);

% Display results
fprintf('Semi-major axis (a): %.4f\n', a);
fprintf('Semi-minor axis (b): %.4f\n', b);

% Plot the ellipse
theta = linspace(0, 2*pi, 100);
x_ellipse = a * cos(theta);
y_ellipse = b * sin(theta);

% Plot ellipse and points
figure;
plot(x_ellipse, y_ellipse, 'b', 'LineWidth', 1.5); % Ellipse
hold on;

% Plot the foci
foci = plot(F1(1), F1(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % F1
plot(F2(1), F2(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % F2

% Plot point P
point_p = plot(Px, Py, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % P

% Configure axis and grid
axis equal;
grid on;
xticks(-6:1:6);
yticks(-4:1:4);

% Add annotations
title(['F = ', num2str(F), ...
       ', P = [', num2str(Px), ',', num2str(Py), ...
       '], a = ', num2str(a), ', b = ', num2str(b)]);
xlabel('x');
ylabel('y');

% Add a legend
legend([point_p, foci], {'P', 'F_1, F_2'}, 'Location', 'best');



%anageo2_9

clear
% Parameters for the ellipse
x0 = 3; % Ellipse center x-coordinate
y0 = 2; % Ellipse center y-coordinate
a = 4;  % Semi-major axis
b = 3;  % Semi-minor axis

% Parameters for the line
k = 1;  % Slope of the line
p = -1; % Intercept of the line

% Substitute y = kx + p into the ellipse equation and simplify
% The ellipse equation becomes:
% ((x - x0)^2 / a^2) + (((k*x + p) - y0)^2 / b^2) - 1 = 0
% Expand to a quadratic form: Ax^2 + Bx + C = 0

A = (1 / a^2) + (k^2 / b^2);
B = (2 * k * (p - y0) / b^2) - (2 * x0 / a^2);
C = ((x0^2 / a^2) + ((p - y0)^2 / b^2) - 1);

% Solve the quadratic equation Ax^2 + Bx + C = 0 using the quadratic formula
discriminant = B^2 - 4 * A * C;

if discriminant < 0
    % No real solutions (line does not intersect the ellipse)
    fprintf('The line does not intersect the ellipse.\n');
    x_vals = [];
    y_vals = [];
else
    % Calculate the two roots of the quadratic equation
    x1 = (-B + sqrt(discriminant)) / (2 * A);
    x2 = (-B - sqrt(discriminant)) / (2 * A);
    x_vals = [x1, x2];
    y_vals = k * x_vals + p; % Corresponding y-values
    fprintf('Intersection points:\n');
    for i = 1:length(x_vals)
        fprintf('(%.4f, %.4f)\n', x_vals(i), y_vals(i));
    end
end

% Plot the ellipse
theta = linspace(0, 2*pi, 200);
x_ellipse = x0 + a * cos(theta);
y_ellipse = y0 + b * sin(theta);

figure;
plot(x_ellipse, y_ellipse, 'b', 'LineWidth', 1.5); % Ellipse
hold on;

% Plot the line
x_line = linspace(min(x_ellipse) - 2, max(x_ellipse) + 2, 100);
y_line = k * x_line + p;
plot(x_line, y_line, 'r', 'LineWidth', 1.5); % Line

% Plot intersection points, if any
if ~isempty(x_vals)
    plot(x_vals, y_vals, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % Intersection points
end

% Formatting
axis equal;
grid on;

% Add title only if there is no intersection
if isempty(x_vals)
    title('Eivät leikkaa');
end

xlabel('x');
ylabel('y');



%sini_1

clear
% Given values
m = 2;         % mass in kg
k = 10;        % spring constant in N/m
x0 = 2;        % initial position in meters
v0 = 3;        % initial velocity in m/s

% Step 1: Calculate angular frequency omega
omega = sqrt(k / m);  % omega = sqrt(k/m)

% Step 2: Calculate the period T
T = 2 * pi / omega;   % period = 2*pi / omega

% Step 3: Calculate the amplitude A
A = sqrt(x0^2 + (v0 / omega)^2);  % amplitude calculation

% Step 4: Calculate the phase angle phi
phi = atan2(x0 * omega, v0);  % phase angle using atan2

% Step 5: Create time array from 0 to 3T
t = linspace(0, 3 * T, 500);

% Step 6: Calculate the position x(t) at each time point
x_t = A * sin(omega * t + phi);

% Step 7: Plot the position vs time
figure;
plot(t, x_t, 'b', 'LineWidth', 1.5);  % Plotting the displacement
title(['m = ', num2str(m), ', k = ', num2str(k), ', x_0 = ', num2str(x0), ', v_0 = ', num2str(v0), ...
       ', \\omega = ', num2str(omega, '%.4f'), ', A = ', num2str(A, '%.4f'), ', \\phi = ', num2str(phi, '%.5f')]);

% Adjust title position (moving it slightly upwards)
title(['m = ', num2str(m), ', k = ', num2str(k), ', x_0 = ', num2str(x0), ', v_0 = ', num2str(v0), ...
       ', \\omega = ', num2str(omega, '%.4f'), ', A = ', num2str(A, '%.4f'), ', \\phi = ', num2str(phi, '%.5f')]);
% Adjust ticks to prevent overlapping
xticks(0:1:8);
yticks(-2.5:0.5:2.5);

% Labels and grid
xlabel('aika t'); % Time (s)
ylabel('paikka x'); % Position (m)
grid on;

% Set axis limits
xlim([0, 8]);
ylim([-2.5, 2.5]);

% Display calculated values
fprintf('Calculated angular frequency (omega): %.4f rad/s\n', omega);
fprintf('Calculated period (T): %.4f s\n', T);
fprintf('Calculated amplitude (A): %.4f m\n', A);
fprintf('Calculated phase angle (phi): %.4f rad\n', phi);



%sini_2

clear
% Parameters
R = 4;
L = 0.05;
C = 0.004;
omega = 40*pi;
I = 5;

% Calculate the amplitude (A)
A = sqrt((R*I)^2 + (omega*L*I - I/(omega*C))^2);

% Calculate the phase (theta)
theta = atan2(omega*L*I - I/(omega*C), R*I);

% Time vector
t = linspace(0, 0.15, 1000);  % Time from 0 to 0.15 seconds, 1000 points

% Voltage across resistor (u_R = RI * sin(omega * t))
u_R = R * I * sin(omega * t);

% Voltage across inductor (u_L = omega * L * I * sin(omega * t + pi/2))
u_L = omega * L * I * sin(omega * t + pi/2);

% Voltage across capacitor (u_C = I / (omega * C) * sin(omega * t - pi/2))
u_C = (I / (omega * C)) * sin(omega * t - pi/2);

% Total voltage (u = u_R + u_L + u_C)
u = u_R + u_L + u_C;

% --- Time-Domain Plot ---
figure;  % Create a new figure for the time-domain plot
plot(t, u_R, 'b', 'DisplayName', 'u_R');
hold on;
plot(t, u_L, 'g', 'DisplayName', 'u_L');
plot(t, u_C, 'r', 'DisplayName', 'u_C');
plot(t, u, 'm', 'DisplayName', 'u');
hold off;

% Labels and Title
xlabel('Aika t (s)');
ylabel('Jännite (V)');
title(sprintf('R = %.1f, L = %.2f, C = %.4f, \\omega = %.1f\\pi, I = %.1f, A = %.4f, \\theta = %.3f', ...
    R, L, C, omega/pi, I, A, theta));
legend('show');
grid on;

% Set axis limits and ticks
xlim([0, 0.15]);
ylim([-35, 35]);
xticks(0:0.05:0.15);
yticks(-30:10:30);

% --- Phasor Diagram ---
% Phasor Components
u_R_magnitude = R * I;                % Voltage magnitude across resistor
u_L_magnitude = omega * L * I;        % Voltage magnitude across inductor
u_C_magnitude = I / (omega * C);      % Voltage magnitude across capacitor
u_LC = u_L_magnitude - u_C_magnitude; % Net reactive voltage (inductor - capacitor)

% Total voltage (phasor sum)
u_x = u_R_magnitude;                  % x-component (real axis, resistive voltage)
u_y = u_LC;                           % y-component (imaginary axis, net reactive voltage)

% Phasor plot
figure;  % Create a new figure for the phasor diagram
hold on;
line([0, u_R_magnitude], [0, 0], 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'u_R');       % Resistor (x-axis)
line([0, 0], [0, u_L_magnitude], 'Color', 'g', 'LineWidth', 2, 'DisplayName', 'u_L');       % Inductor (y-axis up)
line([0, 0], [0, -u_C_magnitude], 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'u_C');      % Capacitor (y-axis down)
line([0, u_x], [0, u_y], 'Color', 'c', 'LineWidth', 2, 'DisplayName', 'u');                 % Total voltage
hold off;

% Formatting
xlabel('Re(V)');
ylabel('Im(V)');
title('Phasor Diagram');
legend('show');
axis equal;
grid on;

% Set axis limits for clarity
xlim([-35, 35]);
ylim([-35, 35]);



%sini_3

clear
% Given parameters
omega = 10 * pi; % Angular frequency
t1 = 0.11; y1 = 3; % First point
t2 = 0.32; y2 = 1; % Second point

% Define a time range for the oscillation
T = 2 * pi / omega; % Period
t = linspace(0, 3 * T, 500); % Time array

% Solve for amplitude (A) and phase (phi) using nonlinear equations
f = @(params) [
    params(1) * sin(omega * t1 + params(2)) - y1; % Equation 1
    params(1) * sin(omega * t2 + params(2)) - y2  % Equation 2
];
initial_guess = [1, 0]; % Initial guesses for [A, phi]
options = optimset('Display', 'off');
solution = fsolve(f, initial_guess, options);

% Extract amplitude and phase
A = solution(1);
phi = solution(2);

% Calculate oscillation curve
y = A * sin(omega * t + phi);

% Plot the oscillation
figure;
plot(t, y, 'b', 'LineWidth', 1.5); % Oscillation curve
hold on;

% Highlight t1, y1 with a red circle
point1 = plot(t1, y1, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % Handle for (t1, y1)

% Highlight t2, y2 with a green circle
point2 = plot(t2, y2, 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g'); % Handle for (t2, y2)

% Set axis limits
xlim([0, 0.6]); % Set x-axis range
ylim([-8, 8]);  % Set y-axis range

% Add grid
grid on; % Major grid
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on'); % Minor grid

% Labels and title
xlabel('aika t'); % Time (s)

% Custom title with exact formatting
title(sprintf('\\omega = 10\\pi, t_1 = %.2f, y_1 = %.0f, t_2 = %.2f, y_2 = %.0f, A = %.4f, \\phi = %.5f', ...
    t1, y1, t2, y2, A, phi));

% Add custom legend using the point handles
legend([point1, point2], {'t_1, y_1', 't_2, y_2'}, 'Location', 'best');



%komp_1

clear
% Define parameters
R1 = 0.4;       % Resistance R1 (Ohms)
L = 0.04;       % Inductance (Henries)
R2 = 1.2;       % Resistance R2 (Ohms)
C = 0.05;       % Capacitance (Farads)
omega = 5 * pi; % Angular frequency (rad/s)
U = 5;          % Voltage magnitude (Volts)

% Calculate impedances
Z1 = R1 + 1i * omega * L;              % Impedance Z1 (resistor + inductor)
Z2 = R2 - 1i / (omega * C);            % Impedance Z2 (resistor + capacitor)
Zs = Z1 + Z2;                          % Series impedance
Zr = (Z1 * Z2) / (Z1 + Z2);            % Parallel impedance

% Calculate currents
Is = U / Zs;                           % Series current
Ir = U / Zr;                           % Parallel current

% Convert impedances and currents to polar form for display
Zs_polar = [abs(Zs), angle(Zs) * 180 / pi]; % |Zs| ∠ θ
Zr_polar = [abs(Zr), angle(Zr) * 180 / pi]; % |Zr| ∠ θ
Is_polar = [abs(Is), angle(Is) * 180 / pi]; % |Is| ∠ θ
Ir_polar = [abs(Ir), angle(Ir) * 180 / pi]; % |Ir| ∠ θ

% Display results
disp(['Zs = ', num2str(Zs_polar(1)), ' ∠ ', num2str(Zs_polar(2)), '°']);
disp(['Zr = ', num2str(Zr_polar(1)), ' ∠ ', num2str(Zr_polar(2)), '°']);
disp(['Is = ', num2str(Is_polar(1)), ' ∠ ', num2str(Is_polar(2)), '°']);
disp(['Ir = ', num2str(Ir_polar(1)), ' ∠ ', num2str(Ir_polar(2)), '°']);

% --- First Plot: Impedances ---
figure;
hold on;

% Phasor lines with DisplayName for legend
h1 = line([0, real(Z1)], [0, imag(Z1)], 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Z_1');
h2 = line([0, real(Z2)], [0, imag(Z2)], 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Z_2');
h3 = line([0, real(Zs)], [0, imag(Zs)], 'Color', 'm', 'LineWidth', 2, 'DisplayName', 'Z_s');
h4 = line([0, real(Zr)], [0, imag(Zr)], 'Color', 'c', 'LineWidth', 2, 'DisplayName', 'Z_r');

% Axes lines
line([-1.5, 2], [0, 0], 'Color', 'k', 'LineWidth', 0.5); % Horizontal axis
line([0, 0], [-1.5, 1.5], 'Color', 'k', 'LineWidth', 0.5); % Vertical axis

hold off;

% Formatting for First Plot
xlabel('Re(Z) (Ohms)');
ylabel('Im(Z) (Ohms)');
title(sprintf('R_1 = %.1f, L = %.2f, R_2 = %.1f, C = %.2f, \\omega = %.1f\\pi, U = %.1f\nZ_s = %.2f ∠ %.1f°, Z_r = %.2f ∠ %.1f°', ...
    R1, L, R2, C, omega / pi, U, Zs_polar(1), Zs_polar(2), Zr_polar(1), Zr_polar(2)));
legend([h1, h2, h3, h4], {'Z_1', 'Z_2', 'Z_s', 'Z_r'}, 'Location', 'northwest'); % Explicit legend
grid on;

% Set axis limits and ticks for a uniform grid
xlim([-1.5, 2]); % Adjust the x-axis limits
ylim([-1.5, 1.5]); % Adjust the y-axis limits
xticks(-1.5:0.5:2); % Set x-ticks at intervals of 0.5
yticks(-1.5:0.5:1.5); % Set y-ticks at intervals of 0.5
axis equal; % Equal scaling for both axes

% --- Second Plot: Voltages and Currents ---
figure;
hold on;

% Phasor lines with DisplayName for legend
h5 = line([0, real(U)], [0, imag(U)], 'Color', 'g', 'LineWidth', 2, 'DisplayName', 'U');
h6 = line([0, real(Is)], [0, imag(Is)], 'Color', 'm', 'LineWidth', 2, 'DisplayName', 'I_s');
h7 = line([0, real(Ir)], [0, imag(Ir)], 'Color', 'c', 'LineWidth', 2, 'DisplayName', 'I_r');

% Axes lines
line([-4, 6], [0, 0], 'Color', 'k', 'LineWidth', 0.5); % Horizontal axis
line([0, 0], [-4, 4], 'Color', 'k', 'LineWidth', 0.5); % Vertical axis

hold off;

% Formatting for Second Plot
xlabel('Re(I, U)');
ylabel('Im(I, U)');
title(sprintf('U = %.1f ∠ %.1f°, I_s = %.2f ∠ %.1f°, I_r = %.2f ∠ %.1f°', ...
    U, 0, Is_polar(1), Is_polar(2), Ir_polar(1), Ir_polar(2)));
legend([h5, h6, h7], {'U', 'I_s', 'I_r'}, 'Location', 'northwest'); % Explicit legend
grid on;
axis equal;

% Set axis limits and ticks for the second plot
xlim([-4, 6]); % Adjust the x-axis limits
ylim([-4, 4]); % Adjust the y-axis limits
xticks(-4:1:6); % Set x-ticks at intervals of 1
yticks(-4:1:4); % Set y-ticks at intervals of 1



%eksjalog_1

clear
% Given values
t1 = 1.2;  % Time t1
U1 = 3.5;  % Voltage U1
t2 = 4.3;  % Time t2
U2 = 1.1;  % Voltage U2

% Calculate tau
tau = (t2 - t1) / (log(U1) - log(U2));

% Calculate U0
U0 = exp(log(U1) + t1 / tau);

% Time vector for plotting
t = linspace(0, 5 * tau, 1000); % Plot over the range 0 to 5*tau

% Voltage function
U = U0 * exp(-t / tau);

% Plotting
figure;
plot(t, U, 'b', 'LineWidth', 1, 'HandleVisibility', 'off');
hold on;

% Highlight the points (t1, U1) and (t2, U2)
plot(t1, U1, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 't_1, U_1');
plot(t2, U2, 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g', 'DisplayName', 't_2, U_2');

% Formatting
xlabel('t');
ylabel('U');
title(sprintf('t_1 = %.1f, U_1 = %.1f, t_2 = %.1f, U_2 = %.1f\nU_0 = %.4f, \\tau = %.4f', ...
    t1, U1, t2, U2, U0, tau), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northwest');
grid on;
hold off;



%eksjalog_2

clear
% Given values
Ty = 25;       % Ambient temperature
t1 = 8;        % Time t1
T1 = 60;       % Temperature T1
t2 = 20;       % Time t2
T2 = 35;       % Temperature T2

% Calculate k
k = (log(T1 - Ty) - log(T2 - Ty)) / (t2 - t1);

% Calculate T0
T0 = Ty + exp(log(T1 - Ty) + k * t1);

% Time vector for plotting
t = linspace(0, 5 / k, 1000); % Plot range: 0 to 5/k

% Temperature function
T = Ty + (T0 - Ty) * exp(-k * t);

% Plotting
figure;

% Plot T(t) curve
plot(t, T, 'b', 'LineWidth', 1, 'HandleVisibility', 'off'); % Exclude from legend
hold on;

% Highlight the points (t1, T1) and (t2, T2)
plot(t1, T1, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 't_1, T_1');
plot(t2, T2, 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g', 'DisplayName', 't_2, T_2');

% Formatting
xlabel('t', 'Interpreter', 'latex');
ylabel('T', 'Interpreter', 'latex');
title(sprintf('T_y = %.1f, t_1 = %.1f, T_1 = %.1f, t_2 = %.1f, T_2 = %.1f\nT_0 = %.4f, k = %.4f', ...
    Ty, t1, T1, t2, T2, T0, k), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northeast');
grid on;
hold off;


