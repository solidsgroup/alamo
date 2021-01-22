
dirs = $(shell find . -maxdepth 1 -name "output*" -type d )
thumbnails = $(addsuffix /thumbnail.png, ${dirs})
zips = $(addsuffix .tar.gz, ${dirs})

thumbnails: ${thumbnails}
	@echo Done

zips: ${zips}
	@echo ${zips}
	@echo Done

%/thumbnail.png: %/metadata
	visit -nowin -cli -s thumbnail.py $< thumbnail.png

%.tar.gz: %
	tar cvzf $@.tmp $<
	mv $@.tmp $@

