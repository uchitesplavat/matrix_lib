CC=gcc
CFLAGS=-Wall -Werror -Wextra -std=c11 #-fsanitize=address -g
VALGRIND=valgrind --leak-check=full --leak-resolution=med -s
OS:= $(shell uname -s)
LINFLAGS=-lcheck -pthread -lrt -lm -lsubunit -s
#TC_FLAGS += $(shell pkg-config --cflags check) $(shell pkg-config --libs check)

all: clean test gcov_report memory

test: s21_matrix.a

ifeq ($(OS), Darwin)
	gcc s21_matrix_tests.c s21_matrix.a  -o test `pkg-config --cflags --libs check` -lcheck -lpthread #-fsanitize=address
	./test
else
	gcc s21_matrix_tests.c s21_matrix.a -o test $(LINFLAGS)
	./test
endif
s21_matrix.o: s21_matrix.c s21_matrix.h
	gcc $(CFLAGS) $(LINFLAGS) -c s21_matrix.c 
s21_matrix.a: s21_matrix.o
	ar r s21_matrix.a  s21_matrix.o
	ranlib s21_matrix.a
gcov_report: 
ifeq ($(OS), Darwin)
	gcc --coverage s21_matrix_tests.c -o s21_matrix_tests s21_matrix.c `pkg-config --cflags --libs check`
	./s21_matrix_tests
	gcov s21_matrix_tests-s21_matrix_tests.gcda
	lcov -t "s21_matrix_tests" -o s21_matrix_tests.info -c -d .
	genhtml -o report s21_matrix_tests.info
	open ./report/index.html
else
	gcc --coverage s21_matrix_tests.c -o s21_matrix_tests s21_matrix.c $(LINFLAGS)
	./s21_matrix_tests
	gcov s21_matrix_tests-s21_matrix_tests.gcda
	lcov -t "s21_matrix_tests" -o s21_matrix_tests.info -c -d .  
	genhtml -o report s21_matrix_tests.info
	#xdg-open ./report/index.html
endif

memory: test
ifeq ($(OS), Darwin)
	CK_FORK=no leaks --atExit -- ./test
	rm -rf s21_matrix.o
else
	#sudo dnf install valgrind
	$(VALGRIND) --leak-check=full --show-leak-kinds=all ./test
	rm -rf s21_matrix.o
endif


clean:
	-rm -rf *.o && rm -rf *.gcno && rm -rf *.exec
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.info && rm -rf *.gcov
	-rm -rf ./test && rm -rf ./gcov_report
	-rm -rf ./report/
