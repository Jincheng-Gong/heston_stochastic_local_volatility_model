rm -r build
cmake -H. -Bbuild && cmake --build build -j4

# read -p "Press 'y' to run the program, any other key to exit: " a

# case $a in
# (Y | y) 
#     echo "Running HestonSLVMC..."
#     ./build/HestonSLVMC ;;
# (N | n)
#     exit ;;
# esac

# echo "EOF"
