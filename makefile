include compiler

clean:
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
	@git commit -am "Commit"
	@git push origin master

pull:
	@git reset --hard HEAD
	@git pull
