CC = gcc
CCFLAGS = -x c++ -lstdc++ -Wall -Wextra  
bindir = /usr/local/bin/
libdir = /usr/local/lib/DEPTH-CLONE-2.7.8/
lib = par/ water_models/
executables = DEPTH complex-depth compute_deviation predict-binding-site generate_sites
INSTALL_H = src/install_dir.h

all:
	echo "install program in current directory"
	echo "#define EXE_DIR string(\""$(PWD)"\")" > $(INSTALL_H)
	echo "#define LIB_DIR string(\""$(PWD)"\")" >> $(INSTALL_H)
	$(CC) $(CCFLAGS) -fopenmp src/main.cpp -o DEPTH
	$(CC) $(CCFLAGS) src/complex-depth.cpp -o complex-depth
	$(CC) $(CCFLAGS) src/compute_deviation.cpp -o compute_deviation
	$(CC) $(CCFLAGS) src/predict-binding-site.cpp -o predict-binding-site
	$(CC) $(CCFLAGS) src/generate_sites.cpp -o generate_sites
	echo "done"

install:
	echo "install program in " $(bindir)
	echo "#define EXE_DIR string(\""$(bindir)"\")" > $(INSTALL_H)
	echo "#define LIB_DIR string(\""$(libdir)"\")" >> $(INSTALL_H)
	$(CC) $(CCFLAGS) -fopenmp src/main.cpp -o DEPTH
	$(CC) $(CCFLAGS) src/complex-depth.cpp -o complex-depth
	$(CC) $(CCFLAGS) src/compute_deviation.cpp -o compute_deviation
	$(CC) $(CCFLAGS) src/predict-binding-site.cpp -o predict-binding-site
	$(CC) $(CCFLAGS) src/generate_sites.cpp -o generate_sites
	rm -rf $(libdir); mkdir $(libdir); cp -r $(lib) $(libdir)
	cp $(executables) $(bindir)
	echo "done"

clean:
	rm $(executables); rm $(INSTALL_H)

remove:
	cd $(bindir); rm $(executables) 
	rm -rf $(libdir)

