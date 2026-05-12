function [contourPts, contourInt] = nyxus_ordered_contour_2d(points, intensities)
% Reproduce Nyxus's regular 2D contour extraction for a single ROI.
%
% Nyxus builds contours on a one-pixel-padded ROI-local image, preserves
% the padded coordinates, then offsets them back by the ROI origin. That is
% why contour coordinates can be one pixel larger than the raw ROI bounds.

if isempty(points)
    contourPts = zeros(0, 2);
    contourInt = zeros(0, 1);
    return;
end

points = double(points);
intensities = double(intensities(:));

minX = min(points(:, 1));
maxX = max(points(:, 1));
minY = min(points(:, 2));
maxY = max(points(:, 2));
width = maxX - minX + 1;
height = maxY - minY + 1;

padded = zeros(height + 2, width + 2);
for i = 1:size(points, 1)
    x = points(i, 1) - minX + 1;
    y = points(i, 2) - minY + 1;
    padded(y + 1, x + 1) = intensities(i) + 1;
end

contours = gather_multicontour(padded, width, height);
if isempty(contours)
    contourPts = zeros(0, 2);
    contourInt = zeros(0, 1);
    return;
end

contourPts = zeros(0, 2);
contourInt = zeros(0, 1);
for i = 1:numel(contours)
    K = contours{i};
    contourPts = [contourPts; K(:, 1:2) + [minX, minY]]; %#ok<AGROW>
    contourInt = [contourInt; K(:, 3)]; %#ok<AGROW>
end
end

function contours = gather_multicontour(P, w, h)
border = zeros(h + 2, w + 2);
inside = false;

for y = 0:(h + 1)
    for x = 0:(w + 1)
        pos = linear_pos(x, y, w);
        bi = get_at(border, pos, w);
        pi = get_at(P, pos, w);

        if bi ~= 0 && ~inside
            inside = true;
        elseif pi ~= 0 && inside
            continue;
        elseif pi == 0 && inside
            inside = false;
        elseif pi ~= 0 && ~inside
            border = set_at(border, pos, w, pi);

            checkLocationNr = 1;
            startPos = pos;
            counter = 0;
            counter2 = 0;
            neighborhood = [
                -1,      7
                -3 - w,  7
                -w - 2,  1
                -1 - w,  1
                 1,      3
                 3 + w,  3
                 w + 2,  5
                 1 + w,  5
            ];

            while true
                checkPosition = pos + neighborhood(checkLocationNr, 1);
                newCheckLocationNr = neighborhood(checkLocationNr, 2);

                if checkPosition < 0 || checkPosition >= numel(P)
                    break;
                end

                pi2 = get_at(P, checkPosition, w);
                if pi2 ~= 0
                    if checkPosition == startPos
                        counter = counter + 1;
                        if newCheckLocationNr == 1 || counter >= 3
                            inside = true;
                            break;
                        end
                    end

                    checkLocationNr = newCheckLocationNr;
                    pos = checkPosition;
                    counter2 = 0;
                    border = set_at(border, checkPosition, w, pi2);
                else
                    checkLocationNr = 1 + mod(checkLocationNr, 8);
                    if counter2 > 8
                        counter2 = 0; %#ok<NASGU>
                        break;
                    else
                        counter2 = counter2 + 1;
                    end
                end
            end
        end
    end
end

C = zeros(0, 3);
for y = 0:(h + 1)
    for x = 0:(w + 1)
        inte = border(y + 1, x + 1);
        if inte == 0
            continue;
        end

        hasNeig = false;
        if x > 0
            hasNeig = hasNeig || border(y + 1, x) ~= 0;
        end
        if x < w - 1
            hasNeig = hasNeig || border(y + 1, x + 2) ~= 0;
        end
        if y > 0
            hasNeig = hasNeig || border(y, x + 1) ~= 0;
        end
        if y < h - 1
            hasNeig = hasNeig || border(y + 2, x + 1) ~= 0;
        end
        if x > 0 && y > 0
            hasNeig = hasNeig || border(y, x) ~= 0;
        end
        if x < w - 1 && y > 0
            hasNeig = hasNeig || border(y, x + 2) ~= 0;
        end
        if x > 0 && y < h - 1
            hasNeig = hasNeig || border(y + 2, x) ~= 0;
        end
        if x < w - 1 && y < h - 1
            hasNeig = hasNeig || border(y + 2, x + 2) ~= 0;
        end

        if hasNeig
            C(end + 1, :) = [x, y, inte - 1]; %#ok<AGROW>
        end
    end
end

if isempty(C)
    contours = {};
    return;
end

Cfix = C;
for i = 1:size(C, 1)
    p = C(i, :);
    if has_location(Cfix, [p(1), p(2) - 1]) && ...
            has_location(Cfix, [p(1), p(2) + 1]) && ...
            has_location(Cfix, [p(1) - 1, p(2)]) && ...
            has_location(Cfix, [p(1) + 1, p(2)])
        idx = find_location(Cfix, p(1:2));
        Cfix(idx, :) = [];
    end
end

U = Cfix;
contours = {};
while ~isempty(U)
    origin = U(1, :);
    [looplen, S] = check_loop(U, origin);
    if looplen ~= 0
        contours{end + 1} = S; %#ok<AGROW>
    end
    U = remove_locations(U, S(:, 1:2));
end
end

function [looplen, S] = check_loop(R, origin)
U = R;
S = origin;
U = remove_location(U, origin(1:2));
P = zeros(0, 3);
looplen = 0;
pxTip = origin;

while ~isempty(U)
    cands = find_cands(U, pxTip);
    if size(cands, 1) > 1
        P(end + 1, :) = pxTip; %#ok<AGROW>
    end

    cands = prune_cands(cands, pxTip);
    if isempty(cands)
        dist2org = pxTip(1:2) - origin(1:2);
        if abs(dist2org(1)) == 1 || abs(dist2org(2)) == 1
            looplen = looplen + 1;
            return;
        elseif isempty(P)
            looplen = 0;
            return;
        else
            pxTip = P(end, :);
            P(end, :) = [];
        end
    elseif size(cands, 1) == 1
        looplen = looplen + 1;
        pxTip = cands(1, :);
        S(end + 1, :) = pxTip; %#ok<AGROW>
        U = remove_location(U, pxTip(1:2));
    else
        looplen = 0;
        return;
    end
end
end

function cands = find_cands(U, tip)
cands = zeros(0, 3);
for i = 1:size(U, 1)
    dx = abs(U(i, 1) - tip(1));
    dy = abs(U(i, 2) - tip(2));
    if dx + dy <= 1
        cands(end + 1, :) = U(i, :); %#ok<AGROW>
    end
end

if ~isempty(cands)
    return;
end

for i = 1:size(U, 1)
    dx = abs(U(i, 1) - tip(1));
    dy = abs(U(i, 2) - tip(2));
    if dx == 1 && dy == 1
        cands(end + 1, :) = U(i, :); %#ok<AGROW>
    end
end
end

function cands = prune_cands(cands, origin)
if size(cands, 1) <= 1
    return;
end

best = cands(1, :);
for i = 2:size(cands, 1)
    if better_step(origin, cands(i, :), best)
        best = cands(i, :);
    end
end
cands = best;
end

function tf = better_step(origin, c1, c2)
pos1 = dial_pos(c1(1:2) - origin(1:2));
pos2 = dial_pos(c2(1:2) - origin(1:2));
tf = pos1 > pos2;
end

function pos = dial_pos(d)
if d(1) > 0
    if d(2) < 0
        pos = 2;
    elseif d(2) > 0
        pos = -1;
    else
        pos = 1;
    end
elseif d(1) < 0
    if d(2) < 0
        pos = 4;
    elseif d(2) > 0
        pos = -3;
    else
        pos = 5;
    end
else
    if d(2) < 0
        pos = 3;
    elseif d(2) > 0
        pos = -2;
    else
        pos = 0;
    end
end
end

function pos = linear_pos(x, y, w)
pos = x + y * (w + 2);
end

function v = get_at(M, pos, w)
x = mod(pos, w + 2);
y = floor(pos / (w + 2));
v = M(y + 1, x + 1);
end

function M = set_at(M, pos, w, value)
x = mod(pos, w + 2);
y = floor(pos / (w + 2));
M(y + 1, x + 1) = value;
end

function tf = has_location(A, xy)
tf = any(A(:, 1) == xy(1) & A(:, 2) == xy(2));
end

function idx = find_location(A, xy)
idx = find(A(:, 1) == xy(1) & A(:, 2) == xy(2), 1);
end

function A = remove_location(A, xy)
idx = find_location(A, xy);
if ~isempty(idx)
    A(idx, :) = [];
end
end

function A = remove_locations(A, xy)
for i = 1:size(xy, 1)
    A = remove_location(A, xy(i, :));
end
end
