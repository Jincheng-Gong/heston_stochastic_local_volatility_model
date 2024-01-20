rm -r build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . --config Debug

# read -p "Press 'y' to run the program, any other key to exit: " a

# case $a in
# (Y | y) 
#     echo "Running HestonSLVMC..."
#     ./build/HestonSLVMC ;;
# (N | n)
#     exit ;;
# esac

echo "EOF"
