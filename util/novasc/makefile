LFLAGS=-lm
CFLAGS=-I.
NOVAS_FILES=novas.o novascon.o nutation.o solsys1.o eph_manager.o readeph0.o

install:novas.a
	cp -rf novas.a ../lib
	cp -rf *.h ../include

novas.a:$(NOVAS_FILES)
	ar rvs $@ $^

%.out:novas.a %.c
	gcc $(CFLAGS) $< $(@:.out=.c) $(LFLAGS) -o $@

%.o:%.c
	gcc -c $^ -o $@

clean:
	rm -rf *.o *.out *~ *.a

