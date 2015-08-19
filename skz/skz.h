int sk_threshold3(int M_int, double s, double Nd, double sk_lims[]);
int sk_threshold6(int M_int, float s_float, float d_float, float sk_lims_float[]);
int compute_sk_mask(float *z,int nx,int ny,int my,int m,float nd,float skmin,float skmax,int *mask);
void decimate_timeseries(float *z,int nx,int ny,int my);
