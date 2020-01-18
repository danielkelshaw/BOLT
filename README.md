# **BOLT:** Lattice Boltzmann Simulator

### **Compilation Instructions:**
From within the ```BOLT``` directory enter the command:

```g++ --std=c++11 src/*```

This will generate an output file: ```a.out```

**Run the Simulator:** ```./a.out```

### **Altering Parameters:**
The majority of user-defined variables are stored in: ```BOLT/include/input.h``` - upon altering any parameter the code needs to be recompiled before running to ensure changes are reflected in the output.

### **Future Plans:**
The next major step for this project is the inclusion of *Immersed Boundary Conditions* in order to allow for simulation of flow around objects within the channel - adding real value to the simulator.
