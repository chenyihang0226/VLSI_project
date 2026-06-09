SRC = main.cpp Graph.cpp Net.cpp Node.cpp solution.cpp evaluate.cpp
OBJ = $(SRC:.cpp=.o)

CXX = g++
CFLAGS = -Wall -Wextra -Werror -std=c++17 -pedantic -fopenmp
CFLAGS += -g

# 自动检测平台：先探测 rm（MSYS2/MinGW 环境），不存在则用 del（Windows cmd）
RMFOUND := $(shell rm --version 2>NUL)
ifdef RMFOUND
RM = rm -f
TARGET = main
else
RM = del /f /q
TARGET = main.exe
endif

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $(OBJ)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(OBJ)