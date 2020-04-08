LDFLAGS_COMMON = -framework GLUT -framework OpenGL -lstdc++
CFLAGS_COMMON = -c -Wall -I./ -O3
COMPILER = g++

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = blob

SOURCES    = blobs.cpp \
						 STVK.cpp \
						 TRIANGLE_MESH.cpp \
						 TRIANGLE.cpp \
						 EXTRAFUNCTIONS.cpp \
						 WALL.cpp \
						 TENSOR3.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
