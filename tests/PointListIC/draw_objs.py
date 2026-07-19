#!/usr/bin/env python3
"""
Interactive tkinter tool to draw closed polygons for the Alamo "pointlist" IC
(src/IC/PointList.H) and export them as an .xyzo file.

Usage:
    python3 draw_objs.py [--domain xmin xmax ymin ymax] [--unit m] [--grid spacing]

Controls:
    Left click        place a vertex (snaps to nearby vertices / grid)
    n                 finish current polygon and start a new one
    v                 toggle void/solid flag on the CURRENT (in-progress) polygon,
                      or on the LAST finished polygon if nothing is in progress
    u / Backspace     undo last placed vertex (or delete last finished polygon)
    Delete            delete the last finished polygon
    Escape            cancel the in-progress polygon
    Ctrl-s / w        export to .xyzo
    Apply domain      re-reads the xmin/xmax/ymin/ymax/unit fields and rescales

Notes on the file format (see src/IC/PointList.H):
    - Each line is "x y z ObjNum". Points sharing the same (integer) ObjNum form
      one polygon; PointList::Parse starts a new polygon whenever ObjNum jumps by
      more than obj_num_threshold (default 0.5) from the previous line.
    - The C++ side auto-closes each polygon (re-appends the first vertex if the
      last one doesn't match), so we do not need to repeat the first vertex on
      export, but every polygon here must still have >= 3 distinct vertices.
    - Inside/outside uses an even-odd ray-crossing test, so winding direction of
      the vertices does not matter.
    - Void (hole) polygons are NOT encoded in the file. They are flagged via a
      separate ParmParse array: phi.ic.pointlist.invert = 0 1 0 ... (one entry
      per polygon, in file order; 1 = void). This script prints that line.
"""
import argparse
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

SNAP_RADIUS_PX = 10
POINT_RADIUS_PX = 4
SOLID_COLOR = "#2b6cb0"   # blue
VOID_COLOR = "#c0392b"    # red
ACTIVE_COLOR = "#333333"
HIGHLIGHT_COLOR = "#f6c343"
CANVAS_MARGIN_PX = 40


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--domain", type=float, nargs=4,
                   default=[-300e-6, 300e-6, -100e-6, 500e-6],
                   metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
                   help="Initial physical domain extent (default matches "
                        "input.sandwich_rocfire_flame_diffuse)")
    p.add_argument("--unit", type=str, default="m",
                   help="Length unit label for coordinates (default: m)")
    p.add_argument("--grid", type=float, default=None,
                   help="Optional initial grid spacing (same units as --unit); "
                        "enables snap-to-grid if given")
    p.add_argument("--out", type=str, default="random_objs.xyzo",
                   help="Default export filename")
    return p.parse_args(argv)


class Polygon:
    def __init__(self, void=False):
        self.pts = []       # list of (x_phys, y_phys), in order, open (not repeating first)
        self.void = void
        self.closed = False
        # canvas item ids for incremental redraw
        self.line_ids = []
        self.point_ids = []
        self.fill_id = None


class DrawApp:
    def __init__(self, root, domain, unit, grid, out_default):
        self.root = root
        self.xmin, self.xmax, self.ymin, self.ymax = domain
        self.unit = unit
        self.grid = grid
        self.snap_to_grid = tk.BooleanVar(value=grid is not None)
        self.out_default = out_default

        self.polygons = []      # finished + in-progress Polygon objects
        self.current = None     # Polygon currently being drawn, or None

        self._build_ui()
        self._bind_events()
        self.root.after(50, self.redraw)

    # ---------------------------------------------------------------- UI ---
    def _build_ui(self):
        self.root.title("PointList IC polygon drawer")

        main = ttk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True)

        # Control panel
        panel = ttk.Frame(main, padding=8)
        panel.pack(side=tk.RIGHT, fill=tk.Y)

        domf = ttk.LabelFrame(panel, text="Domain", padding=6)
        domf.pack(fill=tk.X, pady=4)
        self.xmin_var = tk.StringVar(value=str(self.xmin))
        self.xmax_var = tk.StringVar(value=str(self.xmax))
        self.ymin_var = tk.StringVar(value=str(self.ymin))
        self.ymax_var = tk.StringVar(value=str(self.ymax))
        self.unit_var = tk.StringVar(value=self.unit)
        for label, var in (("xmin", self.xmin_var), ("xmax", self.xmax_var),
                            ("ymin", self.ymin_var), ("ymax", self.ymax_var),
                            ("unit", self.unit_var)):
            row = ttk.Frame(domf)
            row.pack(fill=tk.X, pady=1)
            ttk.Label(row, text=label, width=6).pack(side=tk.LEFT)
            ttk.Entry(row, textvariable=var, width=12).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(domf, text="Apply domain", command=self.apply_domain).pack(fill=tk.X, pady=(4, 0))

        gridf = ttk.LabelFrame(panel, text="Snap", padding=6)
        gridf.pack(fill=tk.X, pady=4)
        ttk.Checkbutton(gridf, text="Snap to grid", variable=self.snap_to_grid).pack(anchor=tk.W)
        row = ttk.Frame(gridf)
        row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text="spacing", width=6).pack(side=tk.LEFT)
        self.grid_var = tk.StringVar(value=str(self.grid) if self.grid else "")
        ttk.Entry(row, textvariable=self.grid_var, width=12).pack(side=tk.LEFT, fill=tk.X, expand=True)

        actf = ttk.LabelFrame(panel, text="Actions", padding=6)
        actf.pack(fill=tk.X, pady=4)
        ttk.Button(actf, text="New polygon (n)", command=self.new_polygon).pack(fill=tk.X, pady=1)
        ttk.Button(actf, text="Toggle void (v)", command=self.toggle_void).pack(fill=tk.X, pady=1)
        ttk.Button(actf, text="Undo point (u)", command=self.undo).pack(fill=tk.X, pady=1)
        ttk.Button(actf, text="Delete last polygon", command=self.delete_last_polygon).pack(fill=tk.X, pady=1)
        ttk.Button(actf, text="Clear all", command=self.clear_all).pack(fill=tk.X, pady=1)
        ttk.Button(actf, text="Export .xyzo (Ctrl-S)", command=self.export).pack(fill=tk.X, pady=(8, 1))

        listf = ttk.LabelFrame(panel, text="Polygons", padding=6)
        listf.pack(fill=tk.BOTH, expand=True, pady=4)
        self.listbox = tk.Listbox(listf, height=12)
        self.listbox.pack(fill=tk.BOTH, expand=True)

        self.status_var = tk.StringVar(value="Click to start a polygon.")
        ttk.Label(panel, textvariable=self.status_var, wraplength=220,
                  foreground="#555").pack(fill=tk.X, pady=(6, 0))

        # Canvas
        self.canvas = tk.Canvas(main, background="white")
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    def _bind_events(self):
        self.canvas.bind("<Configure>", lambda e: self.redraw())
        self.canvas.bind("<Motion>", self.on_motion)
        self.canvas.bind("<Button-1>", self.on_click)
        self.root.bind("<Key-n>", lambda e: self.new_polygon())
        self.root.bind("<Key-v>", lambda e: self.toggle_void())
        self.root.bind("<Key-u>", lambda e: self.undo())
        self.root.bind("<BackSpace>", lambda e: self.undo())
        self.root.bind("<Delete>", lambda e: self.delete_last_polygon())
        self.root.bind("<Escape>", lambda e: self.cancel_current())
        self.root.bind("<Control-s>", lambda e: self.export())
        self.root.bind("<Key-w>", lambda e: self.export())

    # ------------------------------------------------------ coordinates ---
    def _canvas_size(self):
        w = self.canvas.winfo_width() or 800
        h = self.canvas.winfo_height() or 600
        return w, h

    def _transform(self):
        """Return (scale, ox, oy) such that
        px = ox + (x - xmin) * scale
        py = oy - (y - ymin) * scale   (y flipped: screen y grows downward)
        Aspect-preserving, centered within the canvas."""
        w, h = self._canvas_size()
        avail_w = max(w - 2 * CANVAS_MARGIN_PX, 10)
        avail_h = max(h - 2 * CANVAS_MARGIN_PX, 10)
        dx = max(self.xmax - self.xmin, 1e-30)
        dy = max(self.ymax - self.ymin, 1e-30)
        scale = min(avail_w / dx, avail_h / dy)
        draw_w = dx * scale
        draw_h = dy * scale
        ox = CANVAS_MARGIN_PX + (avail_w - draw_w) / 2
        oy = CANVAS_MARGIN_PX + (avail_h - draw_h) / 2 + draw_h
        return scale, ox, oy

    def phys_to_px(self, x, y):
        scale, ox, oy = self._transform()
        return ox + (x - self.xmin) * scale, oy - (y - self.ymin) * scale

    def px_to_phys(self, px, py):
        scale, ox, oy = self._transform()
        x = self.xmin + (px - ox) / scale
        y = self.ymin + (oy - py) / scale
        return x, y

    # ------------------------------------------------------------ snap ---
    def _all_vertices_px(self, exclude_current_last=False):
        """Yield (px, py, phys_x, phys_y) for every placed vertex (finished
        polygons + current in-progress polygon)."""
        for poly in self.polygons:
            for (x, y) in poly.pts:
                px, py = self.phys_to_px(x, y)
                yield px, py, x, y

    def snap(self, px, py):
        """Snap a raw pixel position to (in priority order): the current
        polygon's first vertex (to make closing easy), any other existing
        vertex, then the grid. Returns (px_snapped, py_snapped, is_close_hit)."""
        best = None
        best_d2 = SNAP_RADIUS_PX ** 2

        # Priority: current polygon's own first point (enables closing).
        if self.current is not None and len(self.current.pts) >= 3:
            fx, fy = self.current.pts[0]
            fpx, fpy = self.phys_to_px(fx, fy)
            d2 = (fpx - px) ** 2 + (fpy - py) ** 2
            if d2 <= best_d2:
                return fpx, fpy, True

        for vpx, vpy, vx, vy in self._all_vertices_px():
            d2 = (vpx - px) ** 2 + (vpy - py) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best = (vpx, vpy)
        if best is not None:
            return best[0], best[1], False

        if self.snap_to_grid.get():
            spacing = self._grid_spacing()
            if spacing:
                x, y = self.px_to_phys(px, py)
                gx = round(x / spacing) * spacing
                gy = round(y / spacing) * spacing
                return self.phys_to_px(gx, gy)

        return px, py, False

    def _grid_spacing(self):
        try:
            v = float(self.grid_var.get())
            return v if v > 0 else None
        except ValueError:
            return None

    # --------------------------------------------------------- actions ---
    def new_polygon(self):
        if self.current is not None:
            if len(self.current.pts) < 3:
                messagebox.showwarning("Incomplete polygon",
                                        "Current polygon needs at least 3 vertices "
                                        "before starting a new one; it will be discarded.")
                self.polygons.remove(self.current)
            else:
                self.current.closed = True
        self.current = Polygon()
        self.polygons.append(self.current)
        self.status_var.set("Drawing new polygon. Click first vertex again to close it.")
        self.redraw()

    def toggle_void(self):
        target = self.current if self.current is not None else (
            self.polygons[-1] if self.polygons else None)
        if target is None:
            return
        target.void = not target.void
        self.redraw()

    def undo(self):
        if self.current is not None and self.current.pts:
            self.current.pts.pop()
            if not self.current.pts:
                self.polygons.remove(self.current)
                self.current = None
        elif self.polygons:
            self.delete_last_polygon()
        self.redraw()

    def cancel_current(self):
        if self.current is not None:
            self.polygons.remove(self.current)
            self.current = None
            self.status_var.set("Current polygon canceled.")
            self.redraw()

    def delete_last_polygon(self):
        if self.current is not None:
            self.polygons.remove(self.current)
            self.current = None
        elif self.polygons:
            self.polygons.pop()
        self.redraw()

    def clear_all(self):
        self.polygons = []
        self.current = None
        self.redraw()

    def apply_domain(self):
        try:
            xmin = float(self.xmin_var.get())
            xmax = float(self.xmax_var.get())
            ymin = float(self.ymin_var.get())
            ymax = float(self.ymax_var.get())
        except ValueError:
            messagebox.showerror("Invalid domain", "xmin/xmax/ymin/ymax must be numbers.")
            return
        if xmax <= xmin or ymax <= ymin:
            messagebox.showerror("Invalid domain", "Require xmax > xmin and ymax > ymin.")
            return
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.unit = self.unit_var.get().strip() or "m"
        self.redraw()

    # --------------------------------------------------------- drawing ---
    def on_motion(self, event):
        if self.current is None:
            return
        px, py, _ = self.snap(event.x, event.y)
        self._draw_rubber_band(px, py)

    def on_click(self, event):
        if self.current is None:
            self.current = Polygon()
            self.polygons.append(self.current)

        px, py, is_close_hit = self.snap(event.x, event.y)

        if is_close_hit and len(self.current.pts) >= 3:
            self.current.closed = True
            finished = self.current
            self.current = None
            self.status_var.set(
                f"Closed polygon with {len(finished.pts)} vertices. "
                "Press 'n' or click to start a new one.")
            self.redraw()
            return

        x, y = self.px_to_phys(px, py)
        self.current.pts.append((x, y))
        self.status_var.set(f"Current polygon: {len(self.current.pts)} vertices.")
        self.redraw()

    def _draw_rubber_band(self, px, py):
        self.canvas.delete("rubber")
        if not self.current or not self.current.pts:
            return
        lx, ly = self.phys_to_px(*self.current.pts[-1])
        self.canvas.create_line(lx, ly, px, py, fill=ACTIVE_COLOR, dash=(3, 2), tags="rubber")
        self.canvas.create_oval(px - 3, py - 3, px + 3, py + 3,
                                 outline=HIGHLIGHT_COLOR, width=2, tags="rubber")

    def redraw(self):
        self.canvas.delete("all")
        self._draw_domain()
        if self.snap_to_grid.get():
            self._draw_grid()

        self.listbox.delete(0, tk.END)
        for i, poly in enumerate(self.polygons):
            tag = "void" if poly.void else "solid"
            state = "open" if (poly is self.current) else "closed"
            self.listbox.insert(tk.END, f"{i}: {len(poly.pts)} pts  [{tag}]  ({state})")
            self._draw_polygon(poly)

    def _draw_domain(self):
        x0, y0 = self.phys_to_px(self.xmin, self.ymin)
        x1, y1 = self.phys_to_px(self.xmax, self.ymax)
        self.canvas.create_rectangle(x1, y1, x0, y0, outline="#999", width=1)
        self.canvas.create_text((x0 + x1) / 2, y0 + 14,
                                 text=f"[{self.xmin:g}, {self.xmax:g}] {self.unit}",
                                 fill="#777", font=("TkDefaultFont", 8))
        self.canvas.create_text(x1 + 14, (y0 + y1) / 2,
                                 text=f"[{self.ymin:g}, {self.ymax:g}] {self.unit}",
                                 fill="#777", font=("TkDefaultFont", 8), angle=90, anchor="center")

    def _draw_grid(self):
        spacing = self._grid_spacing()
        if not spacing:
            return
        x0, y0 = self.phys_to_px(self.xmin, self.ymin)
        x1, y1 = self.phys_to_px(self.xmax, self.ymax)
        import math
        nx0 = math.ceil(self.xmin / spacing)
        nx1 = math.floor(self.xmax / spacing)
        ny0 = math.ceil(self.ymin / spacing)
        ny1 = math.floor(self.ymax / spacing)
        for i in range(nx0, nx1 + 1):
            gx = i * spacing
            px, _ = self.phys_to_px(gx, self.ymin)
            self.canvas.create_line(px, y0, px, y1, fill="#eee")
        for j in range(ny0, ny1 + 1):
            gy = j * spacing
            _, py = self.phys_to_px(self.xmin, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#eee")

    def _draw_polygon(self, poly):
        color = VOID_COLOR if poly.void else SOLID_COLOR
        pts_px = [self.phys_to_px(x, y) for (x, y) in poly.pts]

        if poly.closed and len(pts_px) >= 3:
            flat = [c for p in pts_px for c in p]
            self.canvas.create_polygon(*flat, fill=color, outline=color,
                                        stipple="gray25", width=2)
        elif len(pts_px) >= 2:
            for i in range(len(pts_px) - 1):
                self.canvas.create_line(*pts_px[i], *pts_px[i + 1], fill=color, width=2)

        for (px, py) in pts_px:
            r = POINT_RADIUS_PX
            self.canvas.create_oval(px - r, py - r, px + r, py + r, fill=color, outline="")

        if poly.pts:
            fx, fy = self.phys_to_px(*poly.pts[0])
            self.canvas.create_oval(fx - 6, fy - 6, fx + 6, fy + 6, outline=color, width=1)

    # -------------------------------------------------------- export -----
    def export(self):
        finished = []
        for poly in self.polygons:
            if poly is self.current:
                if len(poly.pts) >= 3:
                    poly.closed = True
                    finished.append(poly)
                # else: silently drop an empty/too-short in-progress polygon
            else:
                finished.append(poly)

        finished = [p for p in finished if len(p.pts) >= 3]
        if not finished:
            messagebox.showwarning("Nothing to export",
                                    "No closed polygons with >= 3 vertices yet.")
            return

        path = filedialog.asksaveasfilename(
            defaultextension=".xyzo",
            initialfile=self.out_default,
            filetypes=[("xyzo files", "*.xyzo"), ("All files", "*.*")])
        if not path:
            return

        with open(path, "w") as f:
            for obj_num, poly in enumerate(finished):
                for (x, y) in poly.pts:
                    f.write(f"{x:.8e} {y:.8e} 0.0 {obj_num}\n")

        invert = [1 if p.void else 0 for p in finished]
        lines = [
            f"phi.ic.pointlist.file.name = {path}",
            f"phi.ic.pointlist.file.unit = {self.unit}",
        ]
        if any(invert):
            lines.append("phi.ic.pointlist.invert    = " + " ".join(str(v) for v in invert))

        summary = "\n".join(lines)
        print("\n" + summary + "\n")
        messagebox.showinfo("Exported",
                             f"Wrote {sum(len(p.pts) for p in finished)} points, "
                             f"{len(finished)} polygon(s) to:\n{path}\n\n"
                             "Add to your input file:\n\n" + summary)
        self.status_var.set(f"Exported {len(finished)} polygon(s) to {path}")


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    root = tk.Tk()
    root.geometry("1100x750")
    app = DrawApp(root, args.domain, args.unit, args.grid, args.out)
    root.mainloop()


if __name__ == "__main__":
    main()
