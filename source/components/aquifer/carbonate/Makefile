FCFLAGS = 

default: gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	gfortran -dynamiclib -o carbonate.dylib gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	
mac: gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	gfortran -dynamiclib -o carbonate.dylib gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	
windows: gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	gfortran -shared -o carbonate.dll gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o

linux: FCFLAGS = -fpic
	
linux: gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o
	gfortran -shared $(FCFLAGS) -o carbonate.so gw_rom.o link_model.o mars_rsm_sub.o read_geochem.o read_param.o

gw_rom.o: gw_rom.f90 mars_rsm_sub.o read_geochem.o
	gfortran $(FCFLAGS) -c gw_rom.f90

link_model.o: link_model.f90
	gfortran $(FCFLAGS) -c link_model.f90

mars_rsm_sub.o: mars_rsm_sub.f90
	gfortran $(FCFLAGS) -c mars_rsm_sub.f90

read_geochem.o: read_geochem.f90
	gfortran $(FCFLAGS) -c read_geochem.f90

read_param.o: read_param.f90
	gfortran $(FCFLAGS) -c read_param.f90

clean:
	rm -f *.o *.mod carbonate.dll carbonate.dylib carbonate.so
