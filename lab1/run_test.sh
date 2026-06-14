#!/bin/bash
# 运行 Generated Benchmarks 下 8 个 case 的 43 路划分测试
# 拓扑约束: FPGA Graph/MFS2

set -e

K=43
NUM_THREADS=0  # 0 表示自动检测 CPU 核心数

BENCH_DIR="ICCAD2021-TopoPart-Benchmarks/Generated Benchmarks"
TOPO_FILE="ICCAD2021-TopoPart-Benchmarks/FPGA Graph/MFS2"

# 尝试将常见 MinGW 路径加入 PATH，避免缺少 libgcc/libstdc++/libgomp
for p in "/c/mingw64/bin" "/mingw64/bin" "/c/Program Files/mingw64/bin"; do
    if [ -d "$p" ]; then
        export PATH="$p:$PATH"
        echo "Added to PATH: $p"
        break
    fi
done

# 检查关键 DLL
for dll in libgcc_s_seh-1.dll libstdc++-6.dll libgomp-1.dll; do
    if ! command -v "$dll" >/dev/null 2>&1; then
        echo "WARNING: $dll not found in PATH. main.exe may fail silently."
    fi
done

# 结果输出目录
mkdir -p results

echo "========================================"
echo "k-way partitioning test"
echo "k = $K"
echo "Topology = $TOPO_FILE"
echo "Working directory: $(pwd)"
echo "========================================"

# 先测试 main.exe 本身是否可执行
echo ""
echo ">>> Checking main.exe ..."
if [ ! -x ./main.exe ]; then
    echo "ERROR: ./main.exe not found or not executable."
    exit 1
fi
./main.exe || true

echo ""
echo ">>> Starting benchmark tests ..."

for i in $(seq 1 8); do
    CASE="case$i"
    BENCH_PATH="$BENCH_DIR/$CASE"

    echo ""
    echo "========================================"
    echo ">>> Running $CASE ..."
    echo "Command: ./main.exe \"$BENCH_PATH\" $NUM_THREADS $K --topo \"$TOPO_FILE\""
    echo "========================================"

    set +e
    ./main.exe "$BENCH_PATH" $NUM_THREADS $K --topo "$TOPO_FILE"
    EXIT_CODE=$?
    set -e

    echo "<<< $CASE finished with exit code: $EXIT_CODE"
done

echo ""
echo "========================================"
echo "All cases done. Results saved to ./results/"
echo "========================================"
