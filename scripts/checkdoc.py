#!/usr/bin/env python3
import sys,os
me = os.path.dirname(__file__)
sys.path.insert(0,me+"/../docs/source/")
import Inputs

docs, total = Inputs.scrapeInputs(root=me+"/../src/",writeFiles=False)

if docs < total:
    raise Exception("Documentation is missing from inputs.")

