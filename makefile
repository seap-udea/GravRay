include compiler.in

BRANCH=$(shell bash .getbranch)
PROGRAMS=whattimeisit.exe whereami.exe whereisit.exe whereisthisasteroid.exe whereisinsky.exe wherewillitbe.exe scenarioof.exe throwrays.exe whereonearth.exe wherewrtapex.exe whereisapex.exe

all:$(PROGRAMS)

cleancrap:
	@echo "Cleaning crap..."
	@find . -name "*~" -exec rm -rf {} \;
	@find . -name "#" -exec rm -rf {} \;

cleanexe:
	@echo "Cleaning executable..."
	@rm -rf *.pyc *.out *.exe *.log

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

branch:
	@echo "Branch: $(BRANCH)"

commit:
	@echo "Commiting changes to branch $(BRANCH)..."
	@-git commit -am "Commit" 
	@-git push origin $(BRANCH)

pull:
	@echo "Getting the lattest changes from branch $(BRANCH)..."
	@-git reset --hard HEAD
	@-git pull origin $(BRANCH)

edit:
	emacs -nw *.cpp tests/*.sh tests/*.py makefile *.py README.md

unpack:
	@echo "Unpacking large kernels..."
	@cat util/kernels/de430bsp/* > util/kernels/de430.bsp
	@echo "Done."

#Programs that depend on objects.hpp
launchmany.exe scenarioof.exe wherewillitbe.exe:objects.hpp

