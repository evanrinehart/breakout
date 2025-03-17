RAYLIB_VERSION = master
RAYLIB_PATH = vendor/raylib/src
RAYLIB_URL = https://github.com/raysan5/raylib.git

breakout: raylib.h libraylib.a breakout.c
	gcc -Wall -g -o breakout -I. breakout.c libraylib.a -lm

libraylib.a : raylib
	$(MAKE) -C raylib/src
	cp raylib/src/libraylib.a .

raylib.h : raylib
	cp raylib/src/raylib.h .

raylib :
	git clone --depth=1 --branch=$(RAYLIB_VERSION) $(RAYLIB_URL) raylib

