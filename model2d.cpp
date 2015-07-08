#include "model2d.h"
#include <math.h>

Model2D::Model2D(int _D, int _Ndir)
{

    D = _D;
    Ndir = _Ndir;

    nr_sigmas = 0;

    for (float s = sigma_min; s <= .4*L2; s+=sigma_stp) nr_sigmas++;

    sigmas = new float[nr_sigmas];

    nr_sigmas = 0;

    for (float s = sigma_min; s <= .4*L2; s+=sigma_stp) sigmas[nr_sigmas++] = s;

    L2 = (int) floor(ceil(D/2.0) + .5);  // radius/samplingStep
    LL = 2 * L2 + 1;

    kernels = new float*[Ndir*nr_sigmas];
    for (int i = 0; i < Ndir*nr_sigmas; ++i) {
        kernels[i] = new float[LL*LL];
    }

    kernels_avg = new float[Ndir*nr_sigmas];

    kernels_hat = new float*[Ndir*nr_sigmas];
    for (int i = 0; i < Ndir*nr_sigmas; ++i) {
        kernels_hat[i] = new float[LL*LL];
    }

    kernels_hat_sum_2 = new float[Ndir*nr_sigmas];

    sigma_min = 1.0;
    sigma_stp = 0.5;


    ////////////////////////////////
    // formation
    ////////////////////////////////
    for (int i = 0; i < Ndir*nr_sigmas; i++) {

        int direc_idx = i % Ndir;
        int scale_idx = i / Ndir;

        float sigx = sigmas[scale_idx];
        float sigy = L2;                      // broader than sigmax
        float ang = direc_idx * (3.14 / Ndir);

        float vx =  ((float) cos(ang));
        float vy =  ((float) sin(ang));

        kernels_avg[i] = 0; // average

        for (int j = 0; j < LL*LL; j++) {

            int xx = j % LL;
            int yy = j / LL;

            float currx = (xx - L2) *   vx  + (yy - L2) * vy;
            float curry = (xx - L2) * (-vy) + (yy - L2) * vx;

            kernels[i][j] = (float) exp( -(  ((pow(currx,2)/(2*pow(sigx,2))) + (pow(curry,2)/(2*pow(sigy,2)))   ) )   );
            kernels_avg[i] += kernels[i][j];

        }

        kernels_avg[i] /= (float)(LL*LL);

        kernels_hat_sum_2[i] = 0;

        for (int j = 0; j < LL*LL; j++) {

            kernels_hat[i][j] = kernels[i][j] - kernels_avg[i];
            kernels_hat_sum_2[i] += pow(kernels_hat[i][j], 2);

        }
    }

}

Model2D::~Model2D()
{

    for (int i = 0; i < Ndir*nr_sigmas; i++) {
        delete [] kernels[i];
    }
    delete [] kernels;
    kernels = 0;


    for (int i = 0; i < Ndir*nr_sigmas; i++) {
        delete [] kernels_hat[i];
    }
    delete [] kernels_hat;
    kernels_hat = 0;


    delete [] kernels_avg;
    kernels_avg = 0;


    delete [] kernels_hat_sum_2;
    kernels_hat_sum_2 = 0;

}

Model2D::save(const char * inimg_path)
{

    // save float image
//    Image4DSimple outimg1;
//    outimg1.setData((unsigned char *)outimg, in_sz[0], in_sz[1], in_sz[2], 1, V3D_FLOAT32);
//    callback.saveImage(&outimg1, outimg_file);
//    if(inimg) {delete inimg; inimg =0;}

    //
    for (int i = 0; i < Ndir*nr_sigmas; ++i) {
        // save each kernel as float image

    }

}
