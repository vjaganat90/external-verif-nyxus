function d = nyxus_min_sqdist_ordered_2d(point, cloud)
% Mirror Pixel2::min_sqdist for ordered 2D contour clouds.
%
% Nyxus intentionally uses a logarithmic search over an ordered contour
% instead of an exhaustive minimum. Verifiers use this helper when the C++
% feature code calls Pixel2::min_sqdist directly.

if isempty(cloud)
    d = 0.0;
    return;
end

point = double(point);
cloud = double(cloud);
n = size(cloud, 1);
a = 0;
b = n;
extremD = sqdist(point, cloud(1, :));
extremI = 0;
step = floor((b - a) / log(b - a));

while true
    for i = (a + step):step:(b - 1)
        dist = sqdist(point, cloud(i + 1, :));
        if extremD > dist
            extremD = dist;
            extremI = i;
        end
    end

    if extremI >= step
        stepL = step;
    else
        stepL = extremI;
    end

    if extremI + step < n
        stepR = step;
    else
        stepR = n - extremI;
    end

    a = extremI - stepL;
    b = extremI + stepR;

    if b - a <= 10
        step = 1;
    else
        step = floor((b - a) / log(b - a));
    end

    if b - a <= 2
        break;
    end
end

d = extremD;
end

function d = sqdist(a, b)
dx = b(1) - a(1);
dy = b(2) - a(2);
d = dx * dx + dy * dy;
end
