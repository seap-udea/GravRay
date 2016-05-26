include compiler

PROGRAMS=whattimeisit.exe whereami.exe whereisit.exe whereisthisasteroid.exe whereisinsky.exe wherewillitbe.exe scenario.exe launchmany.exe

all:$(PROGRAMS)

cleanexe:
	@echo "Cleaning executable..."
	@rm -rf *.pyc *.out *.exe

clean:cleanexe
	@echo "Cleaning..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#" -exec rm -rf {} \;
	@rm -rf *.pyc *.out *.exe *.png *.dat

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
scenario.exe wherewillitbe.exe:objects.hpp

