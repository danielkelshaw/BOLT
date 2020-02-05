#ifndef IBMSUPPORT_H
#define IBMSUPPORT_H

class IBMSupportClass {

    friend class ObjectsClass;
    friend class IBMNodeClass;

public:

	IBMSupportClass() {idx = 0; jdx = 0; diracVal = 0.0;};
    IBMSupportClass(int i, int j, double deltaVal);
    ~IBMSupportClass() {};

private:

    int idx;
    int jdx;
    double diracVal;

};

#endif
