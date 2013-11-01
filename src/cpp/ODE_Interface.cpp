class ODE_Interface {
public:
	virtual ~ODE_Interface(){};
	virtual void derivative(double* u, double t, double* v) = 0;
};
