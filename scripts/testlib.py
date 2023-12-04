import numpy, yt, pylab, os, pandas, sys, math

def integrate(x,y):
    ret = 0
    for i in range(len(x)-1): ret = ret + 0.5*(y[i+1]+y[i]) * (x[i+1]-x[i])
    return ret

def validate(path, 
            outdir,
            vars = [],
            start = [0,0,0],
            end = [1,1,1],
            generate_ref_data = False,
            reference=None,
            tolerance=1E-8):
            
    ds = yt.load(path)
    dim = int(ds.domain_dimensions[0] > 1) + int(ds.domain_dimensions[1] > 1) + int(ds.domain_dimensions[2] > 1)

    if not type(start) == list: raise Exception("Invalid start, ",start)
    if not type(end)   == list: raise Exception("Invalid end, ",end)
    
    if len(start) == 2: start.append(0.0)
    if len(end) == 2:   end.append(0.0)
    if dim == 2: 
        start[2] = 0
        end[2] = 0

    prof = ds.ray(start,end)

    new_df = prof.to_dataframe([("gas","x"),("gas","y"),("gas","z"),*vars])
    
    if not reference:
        reference = "reference/reference-{}d.csv".format(dim)

    if generate_ref_data:
        new_df.to_csv(reference)
        print("Reference data generated, now exiting noisily so we don't accidentally pass the test")
        exit(-1)
    
    ref_df = pandas.read_csv(reference)

    all_ok = True
    if not type(len)==list:
        tolerance = [tolerance] * len(vars)
    elif len(tolerance) != len(vars):
        raise Exception("Wrong number of tolerance values")
    for var, tol in zip(vars,tolerance):
        new_x, new_var = [numpy.array(_x) for _x in zip(*sorted(zip(new_df["x"],new_df[var])))]
        ref_x, ref_var = [numpy.array(_x) for _x in zip(*sorted(zip(ref_df["x"],ref_df[var])))]

        pylab.clf()
        pylab.plot(ref_x,ref_var,color='C0',label='ref')
        pylab.plot(new_x,new_var,color='C1',label='new',linestyle='--')
        pylab.savefig(outdir+"/{}.png".format(var))
        
        err = numpy.sqrt(integrate(ref_x, (numpy.interp(ref_x, new_x, new_var) - ref_var)**2))
        mag = numpy.sqrt(integrate(ref_x, (numpy.interp(ref_x, new_x, new_var) + ref_var)**2))

        relerr = err/mag

        print("{} abs error".format(var),err)
        print("{} rel error".format(var),relerr)
        
        if relerr > tol: all_ok = False
        if math.isnan(relerr): all_ok = False

    if not all_ok:
        print("ERROR: Regression test failed")
        raise(Exception("One or more errors out of tolerance or nan"))
