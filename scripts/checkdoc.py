#!/usr/bin/env python3
import sys,os
me = os.path.dirname(__file__)
sys.path.insert(0,me+"/../docs/source/")
import Inputs

Inputs.scrapeInputsSimple(root=me+"/../src/",writeFiles=False)

