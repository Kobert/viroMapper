viroMapper.o: src/viroMapper.c
	gcc src/viroMapper.c src/ref_arg.c src/map_arg.c src/map_alignment.c src/ref_hash.c src/ref_slidingWindow.c src/ref_math.c src/ref_primes.c src/ref_io.c -Wall -O3 -lm -o viroMapper

