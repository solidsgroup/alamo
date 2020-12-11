
thumbnails = $(addsuffix /thumbnail.png, $(shell find . -name "output*" -type d))

thumbnails: ${thumbnails}
	@echo Done

%/thumbnail.png: %
	@echo $@ $<
	@visit -nowin -cli -s thumbnail.py $< thumbnail
