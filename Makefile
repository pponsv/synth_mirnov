.PHONY : all clean doc compile

all: install

bld:
	meson bld

meson_configure: bld
	meson setup --wipe bld

install: meson_configure
	meson install -C bld

clean: 
	rm -rf bld
	rm -rf bin/*.so
