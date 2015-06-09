# For vi/vim use Cntl-v <tab> to insert a tab character in 'expandtab' mode

push:
	rm -rf Mike-Witt.github.io
	git clone https://github.com/Mike-Witt/Mike-Witt.github.io.git
	cp Makefile     Mike-Witt.github.io/Portland-Physics
	cp *.html       Mike-Witt.github.io/Portland-Physics
	cp *.ipynb      Mike-Witt.github.io/Portland-Physics
	cp *.py         Mike-Witt.github.io/Portland-Physics
	cd Mike-Witt.github.io;\
        git add --all;\
        git commit -m "checkin from `hostname`";\
        git push
	rm -rf Mike-Witt.github.io
