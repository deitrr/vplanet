default:
	-gcc -o vplanet *.c -lm

debug:
	-gcc -g -o vplanet_db *.c -lm -Wno-div-by-zero

opt:
	-gcc -o vplanet *.c -lm -O3

profile:
	-gcc -pg -o vplanet *.c -lm

test:
	@echo 'Compiling vplanet...'
	-gcc -o vplanet *.c -lm
	@echo 'Running vtest...'
	@vtest