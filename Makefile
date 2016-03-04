DEBUG = -Wall

test.x : functions.cpp wsRadius.cpp collision.cpp
	g++ $(DEBUG) -o $@ $^
