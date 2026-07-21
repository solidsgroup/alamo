import random
import math

random.seed(4)

# Domain
X_MIN, X_MAX = -300e-6, 300e-6
Y_MIN, Y_MAX = -100e-6, 500e-6
DOMAIN_AREA = (X_MAX - X_MIN) * (Y_MAX - Y_MIN)
DOMAIN_W = X_MAX - X_MIN
DOMAIN_H = Y_MAX - Y_MIN

EPS = 30e-6

# Minimum clearance between shapes and from the domain boundary. The user
# wants ~half the domain covered by many objects, which isn't reachable
# with a full eps (20e-6) gap without cutting the object count way down --
# so the gap is relaxed to half an eps. It's still enough to keep
# neighboring diffuse boundaries (width ~ eps) from fully merging, just not
# fully separated.
GAP_MARGIN = EPS / 2        # 10e-6, min gap between shapes
BOUNDARY_MARGIN = EPS / 2   # 10e-6, min gap from the domain edge

SIDE_MIN, SIDE_MAX = 50e-6, 200e-6

# Number of grid cells (one shape per cell). With GAP_MARGIN = 10e-6, an
# 8x8 grid (64 shapes) reaches ~50% covered area while giving many more
# objects than a coarser grid (a 3x3 grid alone would hit ~50% with only
# 9 shapes, but that's far fewer objects).
NCOLS, NROWS = 3, 4


def polygon_area(pts):
    # Shoelace formula
    n = len(pts)
    s = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0


def random_star_polygon(n_pts, side_target, radius_lo=0.94, radius_hi=1.0):
    """Generate a random, mildly-irregular polygon (angles sorted, radii
    drawn from a narrow range so the shape stays close to convex/isotropic
    and fills its bounding box efficiently), normalized so its bounding-box
    longest side equals side_target, and recentered on the bounding-box
    center."""
    # Evenly-spaced angles with a small jitter (rather than fully random
    # angles) keep the polygon near-isotropic -- points don't bunch up on
    # one side -- so the bounding box stays close to square and the
    # area/bbox fill ratio stays high and predictable.
    base = 2 * math.pi / n_pts
    angles = sorted((i * base + random.uniform(-0.35 * base, 0.35 * base)) % (2 * math.pi)
                     for i in range(n_pts))
    radii = [random.uniform(radius_lo, radius_hi) for _ in range(n_pts)]
    pts = [(r * math.cos(a), r * math.sin(a)) for a, r in zip(angles, radii)]

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    w = max(xs) - min(xs)
    h = max(ys) - min(ys)
    longest = max(w, h)
    scale = side_target / longest
    pts = [(x * scale, y * scale) for x, y in pts]

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    bx = (max(xs) + min(xs)) / 2
    by = (max(ys) + min(ys)) / 2
    pts = [(x - bx, y - by) for x, y in pts]

    bound_r = max(math.hypot(x, y) for x, y in pts)
    area = polygon_area(pts)
    return pts, bound_r, area


def place_shapes():
    ncols, nrows = NCOLS, NROWS
    usable_w = DOMAIN_W - 2 * BOUNDARY_MARGIN
    usable_h = DOMAIN_H - 2 * BOUNDARY_MARGIN
    cell_w = usable_w / ncols
    cell_h = usable_h / nrows
    side = min(cell_w, cell_h) - GAP_MARGIN
    assert SIDE_MIN <= side <= SIDE_MAX, \
        f"grid cell size gives side={side:.2e}, outside [SIDE_MIN, SIDE_MAX]"

    usable_x0 = X_MIN + BOUNDARY_MARGIN
    usable_y0 = Y_MIN + BOUNDARY_MARGIN

    shapes = []  # list of (cx, cy, bound_r, pts, area)
    covered_area = 0.0

    for row in range(nrows):
        for col in range(ncols):
            cell_x0 = usable_x0 + col * cell_w
            cell_y0 = usable_y0 + row * cell_h

            # Randomize size a bit per-shape but keep it within what the
            # cell can hold (with GAP_MARGIN clearance on all sides).
            side_lo = max(SIDE_MIN, side * 0.97)
            side_hi = min(SIDE_MAX, side, cell_w - GAP_MARGIN, cell_h - GAP_MARGIN)
            side_target = random.uniform(side_lo, side_hi)

            n_pts = random.randint(6, 12)
            pts, bound_r, area = random_star_polygon(n_pts, side_target)

            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            w = max(xs) - min(xs)
            h = max(ys) - min(ys)

            # Slack inside the cell after reserving GAP_MARGIN clearance
            # around the shape; jitter the placement within that slack so
            # the grid doesn't look perfectly regular.
            slack_x = max(0.0, cell_w - GAP_MARGIN - w)
            slack_y = max(0.0, cell_h - GAP_MARGIN - h)
            jitter_x = random.uniform(0.0, slack_x)
            jitter_y = random.uniform(0.0, slack_y)

            cx = cell_x0 + GAP_MARGIN / 2 + jitter_x + w / 2
            cy = cell_y0 + GAP_MARGIN / 2 + jitter_y + h / 2

            shapes.append((cx, cy, bound_r, pts, area))
            covered_area += area

    return shapes, covered_area


def write_file(shapes, filename):
    with open(filename, "w") as f:
        for obj_num, (cx, cy, bound_r, pts, area) in enumerate(shapes):
            for (x, y) in pts:
                f.write(f"{cx + x:.8e} {cy + y:.8e} 0.0 {obj_num}\n")


if __name__ == "__main__":
    shapes, covered_area = place_shapes()
    print(f"Placed {len(shapes)} shapes")
    print(f"Domain area: {DOMAIN_AREA:.4e} m^2")
    print(f"Covered (sharp polygon) area: {covered_area:.4e} m^2 "
          f"({100*covered_area/DOMAIN_AREA:.1f}% of domain)")
    for i, (cx, cy, bound_r, pts, area) in enumerate(shapes):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        w = max(xs) - min(xs)
        h = max(ys) - min(ys)
        print(f"  shape {i}: center=({cx:.3e},{cy:.3e}) n_pts={len(pts)} "
              f"bbox=({w:.3e},{h:.3e}) longest_side={max(w,h):.3e} area={area:.3e}")
    write_file(shapes, "random_objs.xyzo")
