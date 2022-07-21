
outputs = $(shell find . -maxdepth 1 -type d -name "output*")
pngs = $(outputs:=/output.png)

default: ${pngs}
	@echo Done

%/output.png: %
	-visit -cli -nowin -s plot.py $<


