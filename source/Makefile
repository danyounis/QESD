FC=gfortran# 9.4.0
PRECISION=DOUBLE

MOD_OBJS=prec.o math.o optimize.o

all: task2 task3 task4 task4a2 task5a2 task5m task5v1 task5v2 task6 task6a2 task7 tidy

task2: $(MOD_OBJS) task2.f08
	$(FC) $^ -o task2 -fopenmp

task3: $(MOD_OBJS) task3.f08
	$(FC) $^ -o task3 -fopenmp

task4: $(MOD_OBJS) task4.f08
	$(FC) $^ -o task4 -fopenmp

task4a2: $(MOD_OBJS) task4a2.f08
	$(FC) $^ -o task4a2 -fopenmp

task5a2: $(MOD_OBJS) task5a2.f08
	$(FC) $^ -o task5a2 -fopenmp

task5m: $(MOD_OBJS) task5m.f08
	$(FC) $^ -o task5m -fopenmp

task5v1: $(MOD_OBJS) task5v1.f08
	$(FC) $^ -o task5v1 -fopenmp

task5v2: $(MOD_OBJS) task5v2.f08
	$(FC) $^ -o task5v2 -fopenmp

task6: $(MOD_OBJS) task6.f08
	$(FC) $^ -o task6 -fopenmp

task6a2: $(MOD_OBJS) task6a2.f08
	$(FC) $^ -o task6a2 -fopenmp

task7: $(MOD_OBJS) task7.f08
	$(FC) $^ -o task7 -fopenmp

prec.o: prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: math.f08
	$(FC) -c $< -fopenmp

optimize.o: optimize.f08
	$(FC) -c $<

trial1_ghz_esd.o: trial1_ghz_esd.f08
	$(FC) -c $< -fopenmp

trial2.o: trial2.f08
	$(FC) -c $< -fopenmp

tidy:
	rm -rf *.o *.mod *.mod0
	@echo Build successful.

clean:
	rm -rf *.o *.mod *.mod0 task2 task3 task4 task4a2 task5a2 task5m task5v1 task5v2 task6 task6a2 task7

wipe:
	rm -rf *.dat
