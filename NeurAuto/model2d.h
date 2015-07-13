#ifndef MODEL2D_H
#define MODEL2D_H

class Model2D
{

public:


    int         D;                      // scale
    int         Ndir;                   // to cover semi-circle
    float *     sigmas;                 // multiscale
    int         nr_sigmas;              // number of scales

    int         L2, LL;                 // ...

    float **    kernels;                // ...
    float *     kernels_avg;            // ...
    float **    kernels_hat;            // ...
    float *     kernels_hat_sum_2;      // ...

    static float sigma_min;
    static float sigma_stp;

    Model2D(int _D, int _Ndir);
    ~Model2D();

    void        save();

};

#endif // MODEL2D_H
