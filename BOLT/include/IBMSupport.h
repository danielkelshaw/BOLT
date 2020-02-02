#ifndef IBMSUPPORT_H
#define IBMSUPPORT_H

class IBMSupportClass {

	friend class ObjectsClass;
	friend class IBMNodeClass;

public:

	IBMSupportClass(int i, int j, double deltaVal);
	~IBMSupportClass() {};

private:

	int idx;
	int jdx;
	double diracVal;

};

#endif
