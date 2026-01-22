mkdir -p build
cmake -S . -B build
cmake --build build -j
./build/clinvar_tool --help || true
