CC = g++
LINK = gfortran
CFLAGS = -c -g -O3 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -Wno-write-strings
LIBS = -lX11 -lcpgplot -lpgplot -lgmp -lpng -lstdc++
TARGET = ffasearch
all: $(TARGET) strip clean

FILES := $(wildcard *.cpp)
OBJS = $(FILES:.cpp=.o)

$(TARGET): $(OBJS)
	$(LINK) $(OBJS) -o $@ $(LIBS)

%.o: %.cpp *.hpp *.h
	$(CC) $(CFLAGS) $< -o $@

strip:
	strip $(TARGET)

clean:
	rm -f $(OBJS) *~
