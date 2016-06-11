include compiler.in

PROGRAMS=whattimeisit.exe whereami.exe whereisit.exe whereisthisasteroid.exe whereisinsky.exe wherewillitbe.exe scenarioof.exe throwrays.exe whereonearth.exe

all:$(PROGRAMS)

cleancrap:
	@echo "Cleaning crap..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#" -exec rm -rf {} \;

cleanexe:
	@echo "Cleaning executable..."
	@rm -rf *.pyc *.out *.exe

clean:cleancrap cleanexe
	@echo "Cleaning..."
	@rm -rf *.png *.dat

%.exe:%.opp
	$(CPP) $^ $(LFLAGS) -o $@

%.opp:%.cpp
	$(CPP) -c $(CFLAGS) $^ -o $@

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@

%.o:%.c
	$(CC) -c $(CFLAGS) $^ -o $@

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin master

pull:
	@-git reset --hard HEAD
	@-git pull

edit:
	emacs -nw *.cpp tests/*.sh tests/*.py makefile *.py README.md

unpack:
	@echo "Unpacking large kernels..."
	@cat util/kernels/de430bsp/* > util/kernels/de430.bsp
	@echo "Done."

#Programs that depend on objects.hpp
launchmany.exe scenarioof.exe wherewillitbe.exe:objects.hpp

