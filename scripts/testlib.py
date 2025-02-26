import numpy, yt, pylab, os, pandas, sys, math

yt.set_log_level(50)  # Disable logging

# Computes the numerical integration of a dataset using the trapezoidal rule.
def integrate(x, y):
    ret = 0
    for i in range(len(x) - 1): 
        ret += 0.5 * (y[i + 1] + y[i]) * (x[i + 1] - x[i])
    return ret

# Reads metadata from the simulation output's Header file, extracting variable names, dimensions, and domain information.
def readHeader(path):
    ret = dict()
    f = open(path + "/Header")
    f.readline().replace('\n', '')  # HyperClaw
    nvars = int(f.readline().replace('\n', ''))  # Number of variables
    ret["vars"] = [f.readline().replace('\n', '') for _ in range(nvars)]
    ret["dim"] = int(f.readline().replace('\n', ''))
    ret["time"] = float(f.readline().replace('\n', ''))
    ret["amrlevels"] = int(f.readline().replace('\n', ''))
    ret["geom_lo"] = [float(x) for x in f.readline().replace('\n', '').split()]
    ret["geom_hi"] = [float(x) for x in f.readline().replace('\n', '').split()]
    return ret

def warning(msg):
    print(msg, file=sys.stderr)

# Main function to validate simulation output against reference data, supporting 1D and 2D extractions.
# This function validates simulation output by extracting 1D or 2D slices from the dataset
# and comparing them to reference data using interpolation and error analysis.
def validate(path,
             outdir,
             plot_type=["TwoDContourPlot"],  # Default to 2D for backward compatibility
             vars=[],
             start=[0, 0, 0],
             end=[1, 1, 1],
             axis=2,
             intercept=0,
             generate_ref_data=False,
             reference=None,
             tolerance=1E-8,
             abs_tolerance=None,
             coord='x',
             tight_layout=True,
             filename_suffix=None):

    info = readHeader(path)

    # Load the dataset using yt. This automatically detects the data format and structure.
    ds = yt.load(path)  # Load dataset

    # Load the dataset using yt. This automatically detects the data format and structure.
    dim = sum(dim_size > 1 for dim_size in ds.domain_dimensions)

    if not isinstance(start, list): raise Exception("Invalid start format:", start)
    if not isinstance(end, list): raise Exception("Invalid end format:", end)

    if len(start) == 2: start.append(0.0)
    if len(end) == 2: end.append(0.0)
    if dim == 2: start[2], end[2] = 0, 0

    if not isinstance(plot_type, list):
        plot_type = [plot_type]

    if not reference:
        reference = "reference/reference-{}d.csv".format(dim)

    if generate_ref_data:
        new_df.to_csv(reference)
        print("Reference data generated, now exiting noisily so we don't accidentally pass the test")
        exit(-1)

    if not isinstance(tolerance, list):
        tolerance = [tolerance] * len(vars)
    elif len(tolerance) != len(vars):
        raise Exception("Wrong number of tolerance values")

    if not abs_tolerance:
        abs_tolerance = tolerance
    elif not isinstance(abs_tolerance, list):
        abs_tolerance = [abs_tolerance] * len(vars)
    elif len(abs_tolerance) != len(vars):
        raise Exception("Wrong number of abs_tolerance values")


    # Loop through each variable to compare extracted simulation data against the reference.
    all_ok = True
    for i, (var, tol, abs_tol) in enumerate(zip(vars, tolerance, abs_tolerance)):
        # Load reference dataset before validation
        if isinstance(reference, list):
            ref_df = pandas.read_csv(reference[i])
        else:
            ref_df = pandas.read_csv(reference)

        ref_df = ref_df.sort_values(by=coord)

        if "OneDLinePlot" in plot_type:
            prof = ds.ray(start, end)
            new_df = prof.to_dataframe([("gas", "x"), *vars])
            ref_df = ref_df[[coord, var]]  # Use only necessary columns for 1D validation
        
        elif "TwoDContourPlot" in plot_type:
            slice2d = ds.slice(axis, intercept)
            prof = ds.ray(start, end)
            new_df = prof.to_dataframe([("gas", "x"), ("gas", "y"), ("gas", "z"), *vars])

            len_x = info["geom_hi"][0] - info["geom_lo"][0]
            len_y = info["geom_hi"][1] - info["geom_lo"][1]

            data2d = slice2d.to_frb(width=len_x, height=len_y, resolution=(1000, 1000 * len_y / len_x))[var]

            pylab.clf()
            im = pylab.imshow(data2d, origin='lower', cmap='jet',
                              extent=[info["geom_lo"][0], info["geom_hi"][0], info["geom_lo"][1], info["geom_hi"][1]])
            pylab.colorbar(im)
            pylab.plot([start[0], end[0]], [start[1], end[1]], linestyle='-', color='white')
            if tight_layout:
                pylab.tight_layout()
            pylab.savefig("{}/2d_{}.png".format(outdir, var))

        new_df = new_df.sort_values(by=coord)

        if "OneDLinePlot" in plot_type:
            new_coord, new_var = new_df["x"], new_df[var]
            ref_coord, ref_var = ref_df["x"], ref_df[var]
        elif "TwoDContourPlot" in plot_type:
            new_coord, ref_coord = new_df[coord], ref_df[coord]
            new_var, ref_var = new_df[var], ref_df[var]

        final_coord = sorted(numpy.concatenate((ref_coord, new_coord), axis=None))
        ref_final_var = numpy.interp(final_coord, ref_coord, ref_var)
        new_final_var = numpy.interp(final_coord, new_coord, new_var)

        pylab.clf()
        pylab.plot(final_coord, ref_final_var, color='C0', label='Reference')
        pylab.plot(final_coord, new_final_var, color='C1', label='New', linestyle='--', marker='o', markerfacecolor='None')
        pylab.xlabel(coord)
        pylab.ylabel(var)
        pylab.legend()
        pylab.grid()
        pylab.savefig(f"{outdir}/{var}.png")

        err = numpy.sqrt(integrate(final_coord, (new_final_var - ref_final_var) ** 2))
        mag = numpy.sqrt(integrate(final_coord, (abs(new_final_var) + abs(ref_final_var)) ** 2))
        relerr = err / mag

        print(f"{var} abs error [tolerance={tolerance[i]}]: {err}")
        print(f"{var} norm: {mag}")
        print(f"{var} rel error [tolerance={tolerance[i]}]: {relerr}")

        if relerr > tolerance[i] or math.isnan(relerr):
            all_ok = False

    if not all_ok:
        raise Exception("Regression test failed: one or more errors out of tolerance.")

