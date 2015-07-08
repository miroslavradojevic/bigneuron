#include "tracer.h"
#include <math.h>
#include <fstream>
#include "toolbox.h"
#include <v3d_interface.h>
#include "v3d_message.h"
#include "basic_surf_objs.h"

float Tracer::gcsstd_min = 1.0;
float Tracer::gcsstd_step = 0.5;
float Tracer::step = 1.5;

Tracer::Tracer(int _Niterations, int _Nstates, int _neuron_radius, bool _is2D, float _zDist, int _Ndirs) {

    Niterations = _Niterations;
    Nstates = _Nstates;
    Ndirs = _Ndirs;

    is2D = _is2D;
    zDist = _zDist;

    // if the trace is in 2d then 1 orthogonal, otherwise two orthogonals
    U2 = _neuron_radius;
    U = 2*_neuron_radius+1; // 1st orthogonal

    if (is2D) {
        W2 = 0;
        W = 1;
    }
    else {
        W2 = _neuron_radius;
        W = 2*_neuron_radius+1;
    }

    V2 = 2; // fixed
    V = 2*V2 + 1;

    // gcsstd define
    int cnt = 0;
    for (float sg = gcsstd_min; sg <= 0.4*U2; sg+=gcsstd_step) cnt++;
    gcsstd_nr = cnt;

    gcsstd = new float[gcsstd_nr];
    cnt = 0;
    for (float sg = gcsstd_min; sg <= 0.4*U2; sg+=gcsstd_step) gcsstd[cnt++] = sg;

    img_vals = new float[U*W*V];

    // tt, tta templates define
    tt = new float*[gcsstd_nr]; // cross-section template
    for (int i = 0; i < gcsstd_nr; ++i) {
        tt[i] = new float[U*W*V];
    }

    tta = new float[gcsstd_nr];

    for (int i = 0; i < gcsstd_nr; ++i) {

        float ag = 0;

        // indexing is differnt in 2d and 3d
        if (is2D) {

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    //for (int ww = -W2; ww <= W2; ++ww) { // W2=0, W=1
                        float value = exp(-(pow(uu,2)  )/(2*pow(gcsstd[i],2))); // +pow(ww,2)
                        int v1 = vv + V2; // 0-V
                        int u1 = uu + U2; // 0-U
//                        int w1 = ww + W2; // 0-W
                        tt[i][v1*U+u1] = value;
                        ag += value;
                    //}
                }
            }

        }
        else {

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    for (int ww = -W2; ww <= W2; ++ww) {
                        float value = exp(-(pow(uu,2)+pow(ww,2))/(2*pow(gcsstd[i],2)));
                        int v1 = vv + V2; // 0-V
                        int u1 = uu + U2; // 0-U
                        int w1 = ww + W2; // 0-W
                        tt[i][v1*U*W+w1*U+u1] = value;
                        ag += value;
                    }
                }
            }

        }

        tta[i] = ag/(U*W*V);

    }

    xt = new float**[Niterations];
    for (int i = 0; i < Niterations; ++i) {
        xt[i] = new float*[Nstates];
        for (int j = 0; j < Nstates; ++j) {
            xt[i][j] = new float[7];
        }
    }

    wt = new float*[Niterations];
    for (int i = 0; i < Niterations; ++i) {
        wt[i] = new float[Nstates];
    }

    prevt = new int*[Niterations];
    for (int i = 0; i < Niterations; ++i) {
        prevt[i] = new int[Nstates];
    }



    xc = new float*[Niterations];
    for (int i = 0; i < Niterations; ++i) {
        xc[i] = new float[3];
    }

    rc = new float[Niterations];

    iter_counter = -1;
    last_queue_element_checked = -1;

    p0       = new float*[Nstates];
    for (int i = 0; i < Nstates; ++i) {
        p0[i] = new float[3];
    }

    v0 = new float*[Nstates];
    for (int i = 0; i < Nstates; ++i) {
        v0[i] = new float[3];
    }

    pred_w0  = new float[Nstates*gcsstd_nr];

    p1       = new float*[Ndirs];
    for (int i = 0; i < Ndirs; ++i) {
        p1[i] = new float[3];
    }

    v1 = new float*[Ndirs];
    for (int i = 0; i < Ndirs; ++i) {
        v1[i] = new float[3];
    }

    pred_w1  = new float[Ndirs*gcsstd_nr];

    trans_xyz0   = new float*[1*Nstates*gcsstd_nr];
    for (int i = 0; i < 1*Nstates*gcsstd_nr; ++i) {
        trans_xyz0[i] = new float[7];
    }

    pties0       = new float[1*Nstates*gcsstd_nr];
    lhoods0      = new float[1*Nstates*gcsstd_nr];
    sort_idx0    = new int[1*Nstates*gcsstd_nr];
    check0       = new bool[1*Nstates*gcsstd_nr];

    trans_xyz1  = new float*[Nstates*Ndirs*gcsstd_nr];
    for (int i = 0; i < Nstates*Ndirs*gcsstd_nr; ++i) {
        trans_xyz1[i] = new float[7];
    }

    pties1      = new float[Nstates*Ndirs*gcsstd_nr];
    lhoods1     = new float[Nstates*Ndirs*gcsstd_nr];
    sort_idx1   = new int[Nstates*Ndirs*gcsstd_nr];
    check1      = new bool[Nstates*Ndirs*gcsstd_nr];

    // directions
    vx_0 = new float[Nstates];
    vy_0 = new float[Nstates];
    vz_0 = new float[Nstates];

    calculate_dirs(vx_0, vy_0, vz_0, is2D, Nstates);

    vx_1 = new float[Ndirs];
    vy_1 = new float[Ndirs];
    vz_1 = new float[Ndirs];

    calculate_dirs(vx_1, vy_1, vz_1, is2D, Ndirs);

    rot = new float*[3];
    rot[0] = new float[3];
    rot[1] = new float[3];
    rot[2] = new float[3];

}

Tracer::~Tracer() {

    delete [] gcsstd; gcsstd = 0;

    delete [] img_vals; img_vals = 0;

    for (int i = 0; i < gcsstd_nr; ++i) {
        delete [] tt[i];
    }
    delete [] tt; tt = 0;

    delete [] tta;

    for (int i = 0; i < Niterations; ++i) {
        for (int j = 0; j < Nstates; ++j) {
            delete [] xt[i][j]; xt[i][j] = 0;
        }
        delete [] xt[i]; xt[i] = 0;
    }
    delete [] xt; xt=0;

    for (int i = 0; i < Niterations; ++i) {
        delete [] wt[i];
    }
    delete [] wt; wt = 0;

    for (int i = 0; i < Niterations; ++i) {
        delete [] prevt[i]; prevt[i] = 0;
    }
    delete [] prevt; prevt = 0;




    for (int i = 0; i < Niterations; ++i) {
        delete [] xc[i];
    }
    delete [] xc; xc = 0;

    delete [] rc; rc = 0;




    for (int i = 0; i < Ndirs; ++i) {
            delete [] p1[i]; p1[i] = 0;
    }
    delete [] p1; p1 = 0;

    for (int i = 0; i < Ndirs; ++i) {
            delete [] v1[i]; v1[i] = 0;
    }
    delete [] v1; v1 = 0;

    delete [] pred_w1; pred_w1 = 0;


    for (int i = 0; i < Nstates; ++i) {
        delete [] p0[i]; p0[i] = 0;
    }
    delete [] p0; p0 = 0;

    for (int i = 0; i < Nstates; ++i) {
        delete [] v0[i]; v0[i] = 0;
    }
    delete [] v0; v0 = 0;

        for (int i = 0; i < Nstates*Ndirs*gcsstd_nr; ++i) {
            delete [] trans_xyz1[i]; trans_xyz1[i] = 0;
        }
        delete [] trans_xyz1; trans_xyz1 = 0;

        delete [] pties1; pties1 = 0;
        delete [] lhoods1; lhoods1 = 0;
        delete [] sort_idx1; sort_idx1 = 0;
        delete [] check1; check1 = 0;

        for (int i = 0; i < 1*Nstates*gcsstd_nr; ++i) {
            delete [] trans_xyz0[i]; trans_xyz0[i] = 0;
        }
        delete [] trans_xyz0; trans_xyz0 = 0;

        delete [] pties0; pties0 = 0;
        delete [] lhoods0; lhoods0 = 0;
        delete [] sort_idx0; sort_idx0 = 0;
        delete [] check0; check0 = 0;

    delete [] vx_0; vx_0 = 0;
    delete [] vy_0; vy_0 = 0;
    delete [] vz_0; vz_0 = 0;

    delete [] vx_1; vx_1 = 0;
    delete [] vy_1; vy_1 = 0;
    delete [] vz_1; vz_1 = 0;

    delete [] rot[0]; rot[0] = 0;
    delete [] rot[1]; rot[1] = 0;
    delete [] rot[2]; rot[2] = 0;
    delete [] rot;  rot = 0;


}

void Tracer::calculate_dirs(float * _vx, float * _vy, float * _vz, bool is2D, int nr_dirs)
{
    // different schemes depending on the dimensionality
    if (is2D) {
        for (int di = 0; di < nr_dirs; ++di) {
            float ang = di * (3.14 / (nr_dirs-1));
            _vx[di] =  cos(ang);
            _vy[di] =  sin(ang);
            _vz[di] =  0;
        }
    }
    else {
        //
        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < nr_dirs; k++) {

            h_k = -1 + 1 * (double)k/(nr_dirs-1); // -1 : 0
            theta_k = acos(h_k);

                    if(k==0 || k==(nr_dirs-1)){

                        phi_k   = 0;
                        phi_k_1 = 0;

                    }
                    else{

                        phi_k = phi_k_1 + 3.6 / ( sqrt(nr_dirs) * sqrt(1-h_k*h_k));
                        phi_k_1 = phi_k;

                    }

                    _vx[k] = (float) (sin(theta_k) * cos(phi_k));
                    _vy[k] = (float) (sin(theta_k) * sin(phi_k));
                    _vz[k] = (float)  cos(theta_k);

        }

    }

}

void Tracer::rotation_matrix(float a1, float a2, float a3, float b1, float b2, float b3, float** R)
{
    // from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    // assume a(a1,a2,a3) and b(b1,b2,b3) are unit vectors

    // v is cross product of (a1, a2, a3) and (b1, b2, b3)
    float v1 = a2*b3 - b2*a3;
    float v2 = -(a1*b3-b1*a3);
    float v3 = a1*b2-b1*a2;

    // cross product is zero - symmetric operations
    float cross_prod_norm_2     = v1*v1+v2*v2+v3*v3;
    float dot_prod              = a1*b1+a2*b2+a3*b3;

    if (cross_prod_norm_2<=0.00001) {
        if (abs(dot_prod-1)<abs(dot_prod+1)) {
            // identity mapping (a and b are aligned)
            R[0][0] = 1;
            R[0][1] = 0;
            R[0][2] = 0;

            R[1][0] = 0;
            R[1][1] = 1;
            R[1][2] = 0;

            R[2][0] = 0;
            R[2][1] = 0;
            R[2][2] = 1;
        }
        else {
            // inversion (a and b are opposite)
            R[0][0] = -1;
            R[0][1] = 0;
            R[0][2] = 0;

            R[1][0] = 0;
            R[1][1] = -1;
            R[1][2] = 0;

            R[2][0] = 0;
            R[2][1] = 0;
            R[2][2] = -1;
        }
    }
    else {
        // cross product is not very small - it's safe to calculate
        float tt = (1-(a1*b1+a2*b2+a3*b3))/(v1*v1+v2*v2+v3*v3);
        R[0][0] = 1 + 0     + tt * (-v3*v3-v2*v2);
        R[0][1] = 0 + (-v3) + tt * (v1*v2);
        R[0][2] = 0 + (v2)  + tt * (-v1*v3);

        R[1][0] = 0 + (v3)  + tt * (v1*v2);
        R[1][1] = 1 + 0     + tt * (-v3*v3-v1*v1);
        R[1][2] = 0 + (-v1) + tt * (v2*v3);

        R[2][0] = 0 + (-v2) + tt * (v1*v3);
        R[2][1] = 0 + (v1)  + tt * (v2*v3);
        R[2][2] = 1 + 0     + tt * (-v2*v2-v1*v1);
    }

}

void Tracer::rotation_apply(float** R, float v1, float v2, float v3, float &out1, float &out2, float &out3)
{

    out1 = R[0][0]*v1 + R[0][1]*v2 + R[0][2]*v3;
    out2 = R[1][0]*v1 + R[1][1]*v2 + R[1][2]*v3;
    out3 = R[2][0]*v1 + R[2][1]*v2 + R[2][2]*v3;

}

float Tracer::interp(float atX, float atY, float atZ, unsigned char * img, int width, int height, int length)
{

    int x1 = (int) atX;
    int x2 = x1 + 1;
    float x_frac = atX - x1;

    int y1 = (int) atY;
    int y2 = y1 + 1;
    float y_frac = atY - y1;

    if (length==1) { // atZ is not necessary

        bool isIn2D = x1>=0 && x2<width && y1>=0 && y2<height;

        if(!isIn2D) {
            printf("interp() out of boundary [%6.2f, %6.2f, %6.2f],[%d--%d],[%d--%d] M=%d, N=%d, P=%d \n", atX, atY, atZ, x1, x2, y1, y2, width, height, length);
            return 0;
        }

        // take neighbourhood 2d
        float I11_1 = img[y1*width+x1]; // img3d_xyz[ x1 ][ y1 ][ z1 ];  // upper left
        float I12_1 = img[y1*width+x2]; // img3d_xyz[ x2 ][ y1 ][ z1 ];  // upper right
        float I21_1 = img[y2*width+x1]; // img3d_xyz[ x1 ][ y2 ][ z1 ];  // bottom left
        float I22_1 = img[y2*width+x2]; // img3d_xyz[ x2 ][ y2 ][ z1 ];  // bottom right

        return (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1);

    }
    else {

        int z1 = (int) atZ;
        int z2 = z1 + 1;
        float z_frac = atZ - z1;

        bool isIn3D = y1>=0 && y2<height && x1>=0 && x2<width && z1>=0 && z2<length;

        if(!isIn3D) {
            printf("interp() out of boundary [%6.2f, %6.2f, %6.2f],[%d--%d],[%d--%d],[%d--%d] M=%d, N=%d, P=%d \n", atX, atY, atZ, x1, x2, y1, y2, z1, z2, width, height, length);
            return 0;
        }

        // take neighbourhood 3d
        float I11_1 = img[z1*width*height+y1*width+x1]; // img3d_xyz[ x1 ][ y1 ][ z1 ];  // upper left
        float I12_1 = img[z1*width*height+y1*width+x2]; // img3d_xyz[ x2 ][ y1 ][ z1 ];  // upper right
        float I21_1 = img[z1*width*height+y2*width+x1]; // img3d_xyz[ x1 ][ y2 ][ z1 ];  // bottom left
        float I22_1 = img[z1*width*height+y2*width+x2]; // img3d_xyz[ x2 ][ y2 ][ z1 ];  // bottom right

        float I11_2 = img[z2*width*height+y1*width+x1]; // img3d_xyz[ x1 ][ y1 ][ z2 ]; // upper left
        float I12_2 = img[z2*width*height+y1*width+x2]; // img3d_xyz[ x2 ][ y1 ][ z2 ]; // upper right
        float I21_2 = img[z2*width*height+y2*width+x1]; // img3d_xyz[ x1 ][ y2 ][ z2 ]; // bottom left
        float I22_2 = img[z2*width*height+y2*width+x2]; // img3d_xyz[ x2 ][ y2 ][ z2 ]; // bottom right

        return (1-z_frac)  *
                (  (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1) )   +
                        z_frac      *
                (  (1-y_frac) * ((1-x_frac)*I11_2 + x_frac*I12_2) + (y_frac) * ((1-x_frac)*I21_2 + x_frac*I22_2) );

    }

}

void Tracer::trace(float x,  float y,  float z,
                   float vx, float vy, float vz,
                   float gcsstd,
                   unsigned char * img,
                   int  img_width,
                   int img_height,
                   int img_length,
                   float angula_diff_std_deg,
                   float gcsstd_diff_std_pix) {

    // set the start values
    start_px = x;
    start_py = y;
    start_pz = z;
    start_vx = vx;
    start_vy = vy;
    start_vz = vz;

    iter_counter = 0;

    int swcdbg_cnt = 1;

    while (iter_counter<Niterations) { // Niterations

        if (iter_counter == 0) {
            iter_beg(x, y, z, vx, vy, vz, gcsstd, img, img_width, img_height, img_length, angula_diff_std_deg, gcsstd_diff_std_pix);
            for (int i = 0; i < Nstates*gcsstd_nr; ++i) {
                printf("%d\t 6\t %f\t %f\t %f\t %f\t -1 \n",
                       swcdbg_cnt,
                        trans_xyz0[sort_idx0[i]][0],
                        trans_xyz0[sort_idx0[i]][1],
                        trans_xyz0[sort_idx0[i]][2],
                        0.25); // 0.5*pties0[i]/Nstates
                swcdbg_cnt++;
            }
        }
        else {
            iter(iter_counter, img, img_width, img_height, img_length, angula_diff_std_deg, gcsstd_diff_std_pix);
            for (int i = 0; i < Nstates*Ndirs*gcsstd_nr; ++i) {
                printf("%d\t 6\t %f\t %f\t %f\t %f\t -1 \n",
                       swcdbg_cnt,
                        trans_xyz1[sort_idx1[i]][0],
                        trans_xyz1[sort_idx1[i]][1],
                        trans_xyz1[sort_idx1[i]][2],
                        0.25); // 0.5*pties1[i]
                swcdbg_cnt++;
            }
        }

        for (int i = 0; i < Nstates; ++i) {
            printf("%d\t 5\t %f\t %f\t %f\t %f\t -1 \n",
                   swcdbg_cnt,
                    xt[iter_counter][i][0],
                    xt[iter_counter][i][1],
                    xt[iter_counter][i][2],
                    0.5); // wt[iter_counter][i]
            swcdbg_cnt++;
        }

        iter_counter++; // will be used to describe how far we've gone

    }

}

void Tracer::predict1(  float px, float py, float pz,
                        float vx, float vy, float vz,
                        float vx_prev, float vy_prev, float vz_prev,
                        float gcsstd_pv,
                        float angula_diff_std_deg, float gcsstd_diff_std_pix)
{

    for (int i = 0; i < Ndirs; ++i) {

        // set of 2d/3d directions distributed along x/z direction is rotated
        // vx_1,vy_1,vz_1 -> v1[][0,1,2]

        if (is2D) {

            // (0,1,0)->T->(vx,vy,vz)
            rotation_matrix(0,1,0,          vx,vy,vz,                       rot);
            rotation_apply(rot,     vx_1[i],vy_1[i],vz_1[i],    v1[i][0],v1[i][1],v1[i][2]);

            v1[i][2] = 0; // just in case, this is 2d

            p1[i][0] = px + step * v1[i][0]; // px
            p1[i][1] = py + step * v1[i][1]; // py
            p1[i][2] = pz;

            float dot_prod = v1[i][0]*vx+v1[i][1]*vy+v1[i][2]*vz;
            dot_prod = (dot_prod>1)? 1 : (dot_prod<-1)? -1 : dot_prod ;
            float prior_weight = exp(-pow(acos(dot_prod)*(180.0/3.14),2)/(2*pow(angula_diff_std_deg,2)));

            float dvx1 = v1[i][0]-vx;
            float dvy1 = v1[i][1]-vy;
            float dvz1 = v1[i][2]-vz;

            float dvx2 = vx-vx_prev;
            float dvy2 = vy-vy_prev;
            float dvz2 = vz-vz_prev;

            dot_prod = dvx1*dvx2+dvy1*dvy2+dvz1*dvz2;
            dot_prod = (dot_prod>1)? 1 : (dot_prod<-1)? -1 : dot_prod ;
//            prior_weight *= exp(-pow(acos(dot_prod)*(180.0/3.14),2)/(2*pow(angula_diff_std_deg,2)));

            for (int j = 0; j < gcsstd_nr; j++) {
                pred_w1[i*gcsstd_nr+j] =
                        prior_weight *1;// exp(-pow(gcsstd[j]-gcsstd_pv,2)/(2*pow(gcsstd_diff_std_pix,2)));
            }

        }
        else {

            // (0,0,-1)->T->(vx,vy,vz)
            rotation_matrix(0,0,-1,          vx,vy,vz,                       rot);
            rotation_apply(rot,     vx_1[i],vy_1[i],vz_1[i],    v1[i][0],v1[i][1],v1[i][2]);

            p1[i][0] = px +                 step * v1[i][0]; // px
            p1[i][1] = py +                 step * v1[i][1]; // py
            p1[i][2] = pz + (1.0/zDist) *   step * v1[i][2]; // pz

            float dot_prod = v1[i][0]*vx+v1[i][1]*vy+v1[i][2]*vz;
            dot_prod = (dot_prod>1)? 1 : (dot_prod<-1)? -1 : dot_prod ;
            float ang_div_rad = acos(dot_prod);
            float ang_div_deg = ang_div_rad * (180.0/3.14);

            for (int j = 0; j < gcsstd_nr; j++) {
                pred_w1[i*gcsstd_nr+j] =
                        exp(-pow(ang_div_deg,2)/(2*pow(angula_diff_std_deg,2))) *
                        1;//exp(-pow(gcsstd[j]-gcsstd_pv,2)/(2*pow(gcsstd_diff_std_pix,2)));
            }

        }

    } // i .. Ndirs

}

void Tracer::predict0(  float px, float py, float pz,
                        float vx, float vy, float vz,
                        float gcsstd_pv,
                        float angula_diff_std_deg, float gcsstd_diff_std_pix)
{

    for (int i = 0; i < Nstates; ++i) { // this is the first prediction

        // set of 2d/3d directions distributed along x/z direction is rotated
        // vx_0,vy_0,vz_0 -> v0[][0,1,2]

        if (is2D) {

            // (0,1,0)->T->(vx,vy,vz)
            rotation_matrix(0,1,0,          vx,vy,vz,                       rot);
            rotation_apply(rot,     vx_0[i],vy_0[i],vz_0[i],    v0[i][0],v0[i][1],v0[i][2]);

            v0[i][2] = 0; // just in case, this is 2d

            p0[i][0] = px + step * v0[i][0]; // px
            p0[i][1] = py + step * v0[i][1]; // py
            p0[i][2] = pz;

            float dot_prod = v0[i][0]*vx+v0[i][1]*vy+v0[i][2]*vz;
            dot_prod = (dot_prod>1)? 1 : (dot_prod<-1)? -1 : dot_prod ;
            float ang_div_deg = acos(dot_prod) * (180.0/3.14);

            for (int j = 0; j < gcsstd_nr; j++) {
                pred_w0[i*gcsstd_nr+j] =
                        exp(-pow(ang_div_deg,2)/(2*pow(angula_diff_std_deg,2))) *
                        exp(-pow(gcsstd[j]-gcsstd_pv,2)/(2*pow(gcsstd_diff_std_pix,2)));// will affect lhood and pties and eventually centroid calculation
            }

        }
        else {

            // (0,0,-1)->T->(vx,vy,vz)
            rotation_matrix(0,0,-1,          vx,vy,vz,                       rot);
            rotation_apply(rot,     vx_0[i],vy_0[i],vz_0[i],    v0[i][0],v0[i][1],v0[i][2]);

            p0[i][0] = px +                 step * v0[i][0]; // px
            p0[i][1] = py +                 step * v0[i][1]; // py
            p0[i][2] = pz + (1.0/zDist) *   step * v0[i][2]; // pz

            float dot_prod = v0[i][0]*vx+v0[i][1]*vy+v0[i][2]*vz;
            dot_prod = (dot_prod>1)? 1 : (dot_prod<-1)? -1 : dot_prod ;
            float ang_div_deg = acos(dot_prod) * (180.0/3.14);

            for (int j = 0; j < gcsstd_nr; j++) {
                pred_w0[i*gcsstd_nr+j] =
                        exp(-pow(ang_div_deg,2)/(2*pow(angula_diff_std_deg,2))) *
                        exp(-pow(gcsstd[j]-gcsstd_pv,2)/(2*pow(gcsstd_diff_std_pix,2)));
            }

        }

    } // i .. Nstates

}

void Tracer::iter_beg(  float               xloc,
                    float                   yloc,
                    float                   zloc,
                    float                   vxloc,
                    float                   vyloc,
                    float                   vzloc,
                    float                   gcsstdcloc,
                    unsigned char*          img,
                    int                     img_width,
                    int                     img_height,
                    int                     img_length,
                    float                   angula_diff_std_deg,
                    float                   gcsstd_diff_std_pix) {

    predict0(xloc, yloc, zloc, vxloc, vyloc, vzloc, gcsstdcloc, angula_diff_std_deg, gcsstd_diff_std_pix); // p0, v0, pred_w0

    int count = 0;
    float sum_lhoods = 0;// used to check validity
    float sum_prior = 0;

    for (int j = 0; j < Nstates*gcsstd_nr; ++j) {    // calculate likelihoods

        trans_xyz0[count][0] = p0[j/gcsstd_nr][0];
        trans_xyz0[count][1] = p0[j/gcsstd_nr][1];
        trans_xyz0[count][2] = p0[j/gcsstd_nr][2];

        trans_xyz0[count][3] = v0[j/gcsstd_nr][0];
        trans_xyz0[count][4] = v0[j/gcsstd_nr][1];
        trans_xyz0[count][5] = v0[j/gcsstd_nr][2];

        trans_xyz0[count][6] = gcsstd[j%gcsstd_nr];

        lhoods0[count] = interp(trans_xyz0[count][0], trans_xyz0[count][1], trans_xyz0[count][2], img, img_width, img_height, img_length);
//                zncc(
//                        trans_xyz0[count][0], trans_xyz0[count][1], trans_xyz0[count][2],
//                        trans_xyz0[count][3], trans_xyz0[count][4], trans_xyz0[count][5], j%gcsstd_nr,
//                        img, img_width, img_height, img_length);

        sum_lhoods += lhoods0[count];

        pties0[count] = pred_w0[j]; // stores priors

        sum_prior += pties0[count];

        count++;

    } // loop Nstates*gcsstd_nr

    for (int i = 0; i < Nstates*gcsstd_nr; i++) {
        pties0[i] =
                1 *
//                pties0[i] * //(pties0[i]/sum_prior) *
                ((sum_lhoods>0.0001)?lhoods0[i]:1);
    }

    // take best Ns (we'll always have enough to take)
    descending(pties0, Nstates*gcsstd_nr, sort_idx0); // ptes0 will be sorted as a side effect

    // top Nstates will be selected after sorting
    bool check[Nstates*gcsstd_nr]; // static alloc TODO
    for (int var = 0; var < Nstates*gcsstd_nr; ++var) {
        check[var] = false;
    }

    int cnt = 0;
    for (int i = 0; i < Nstates*gcsstd_nr; i++) {

        int curr_max_idx = sort_idx0[i];
        int loc_idx = curr_max_idx/gcsstd_nr; // y coord

        if (check[curr_max_idx]) continue;

        for (int var = 0; var < gcsstd_nr; ++var) {
            check[loc_idx*gcsstd_nr+var] = true; // fill it for all scales at that location
        }

        // 0th iteration there were Nstates predictions
        xt[0][cnt][0] = trans_xyz0[sort_idx0[i]][0];
        xt[0][cnt][1] = trans_xyz0[sort_idx0[i]][1];
        xt[0][cnt][2] = trans_xyz0[sort_idx0[i]][2];
        xt[0][cnt][3] = trans_xyz0[sort_idx0[i]][3];
        xt[0][cnt][4] = trans_xyz0[sort_idx0[i]][4];
        xt[0][cnt][5] = trans_xyz0[sort_idx0[i]][5];
        xt[0][cnt][6] = trans_xyz0[sort_idx0[i]][6];

        wt[0][cnt]    = pties0[i];       // because they are already sorted

        prevt[0][cnt] = -1; // means it's pointing to start_*

        cnt++;

        if (cnt==Nstates) break;         // add top Ns points that are in the circle for further tracing

    }

    probability_distribution(wt[0], Nstates);       // wt will be normalized

    // estimation centroid xc, rc at the end of the iteration
    xc[0][0] = 0;
    xc[0][1] = 0;
    xc[0][2] = 0;
    rc[0]    = 0;

    for (int j = 0; j < Nstates; ++j) {
        xc[0][0]    += wt[0][j] * xt[0][j][0];
        xc[0][1]    += wt[0][j] * xt[0][j][1];
        xc[0][2]    += wt[0][j] * xt[0][j][2];
        rc[0]       += wt[0][j] * xt[0][j][6];
    }



}

void Tracer::iter(  int                     iter,
                    unsigned char *         img,
                    int                     img_width,
                    int                     img_height,
                    int                     img_length,
                    float                   angula_diff_std_deg,
                    float                   gcsstd_diff_std_pix) {

            int count = 0;
            float sum_lhoods = 0;// used to check validity
            float sum_prior = 0;

            for (int i = 0; i < Nstates; i++) { // will loop throughout the states of the previous iteration

                float px        = xt[iter-1][i][0];
                float py        = xt[iter-1][i][1];
                float pz        = xt[iter-1][i][2];

                float vx        = xt[iter-1][i][3];
                float vy        = xt[iter-1][i][4];
                float vz        = xt[iter-1][i][5];

                float gcsstdc   = xt[iter-1][i][6];

                int idx_prev_state = prevt[iter-1][i];

                float pred_vx, pred_vy, pred_vz;

                if (idx_prev_state==-1) {
                    pred_vx = start_vx;
                    pred_vy = start_vy;
                    pred_vz = start_vz;
                }
                else {
                    pred_vx = xt[iter-2][  idx_prev_state  ][3];
                    pred_vy = xt[iter-2][  idx_prev_state  ][4];
                    pred_vz = xt[iter-2][  idx_prev_state  ][5];
                }



                // calling this one means that it's not the first one - there was a previous iteration
                predict1(
                            px, py, pz,
                            vx, vy, vz,
                            pred_vx, pred_vy, pred_vz,
                            gcsstdc, angula_diff_std_deg, gcsstd_diff_std_pix); // p1, v1, pred_w1

                for (int j = 0; j < Ndirs*gcsstd_nr; j++) { // loop all options, calculate likelihoods

                    trans_xyz1[count][0] = p1[j/gcsstd_nr][0]; // px
                    trans_xyz1[count][1] = p1[j/gcsstd_nr][1]; // py
                    trans_xyz1[count][2] = p1[j/gcsstd_nr][2]; // pz

                    trans_xyz1[count][3] = v1[j/gcsstd_nr][0]; // vx
                    trans_xyz1[count][4] = v1[j/gcsstd_nr][1]; // vy
                    trans_xyz1[count][5] = v1[j/gcsstd_nr][2]; // vz

                    trans_xyz1[count][6] = gcsstd[j%gcsstd_nr]; // gcsstd[]

                    lhoods1[count] = interp(trans_xyz1[count][0], trans_xyz1[count][1], trans_xyz1[count][2], img, img_width, img_height, img_length);

//                            zncc(
//                                    trans_xyz1[count][0], trans_xyz1[count][1], trans_xyz1[count][2],
//                                    trans_xyz1[count][3], trans_xyz1[count][4], trans_xyz1[count][5], j%gcsstd_nr,
//                                    img, img_width, img_height, img_length);

                    sum_lhoods += lhoods1[count];

                    pties1[count] = pred_w1[j]; // w stores priors

                    sum_prior += pties1[count];

                    count++;

                } // Ndirs*gcsstd_nr

            } // Nstates

            for (int i = 0; i < Nstates*Ndirs*gcsstd_nr; i++) {
                pties1[i] =
                        wt[iter-1][i/(Ndirs*gcsstd_nr)] *
//                        pties1[i] * //(pties1[i]/sum_prior) *
                        ((sum_lhoods>0.0001)?lhoods1[i]:1);
            }

            // take best Ns (we'll always have enough to take)
            descending(pties1, Nstates*Ndirs*gcsstd_nr, sort_idx1); // pties1 will be sorted as a side effect

            // top Nstates will be selected after sorting
            for (int var = 0; var < Nstates*Ndirs*gcsstd_nr; ++var) check1[var] = false;

            int cnt = 0;
            for (int i = 0; i < Nstates*Ndirs*gcsstd_nr; i++) {

                int curr_max_idx = sort_idx1[i];
                int loc_idx = curr_max_idx/(gcsstd_nr); // 0 .. Nstates*Ndirs

                if (check1[curr_max_idx]) continue;

                for (int var = 0; var < gcsstd_nr; ++var) {
                    check1[loc_idx*gcsstd_nr+var] = true; // fill it for all scales at that location
                }

                xt[iter][cnt][0] = trans_xyz1[sort_idx1[i]][0];
                xt[iter][cnt][1] = trans_xyz1[sort_idx1[i]][1];
                xt[iter][cnt][2] = trans_xyz1[sort_idx1[i]][2];
                xt[iter][cnt][3] = trans_xyz1[sort_idx1[i]][3];
                xt[iter][cnt][4] = trans_xyz1[sort_idx1[i]][4];
                xt[iter][cnt][5] = trans_xyz1[sort_idx1[i]][5];
                xt[iter][cnt][6] = trans_xyz1[sort_idx1[i]][6];

                wt[iter][cnt]    = pties1[i]; // because they are already sorted

                prevt[iter][cnt] = curr_max_idx/(Ndirs*gcsstd_nr); // 0 .. Nstates

                cnt++;

                if (cnt==Nstates) break; // add top Ns points that are in the circle for further tracing

            }

            probability_distribution(wt[iter], Nstates); // _wt will be normalized

            // set xc, rc at the end of the iteration
            xc[iter][0] = 0;
            xc[iter][1] = 0;
            xc[iter][2] = 0;
            rc[iter]    = 0;

            for (int j = 0; j < Nstates; ++j) {
                xc[iter][0]    += wt[iter][j] * xt[iter][j][0];
                xc[iter][1]    += wt[iter][j] * xt[iter][j][1];
                xc[iter][2]    += wt[iter][j] * xt[iter][j][2];
                rc[iter]       += wt[iter][j] * xt[iter][j][6];
            }

}

float Tracer::zncc(float x, float y, float z, float vx, float vy, float vz, int gcsstd_idx, unsigned char * img, int img_w, int img_h, int img_l) {

    // loop as it was 1:1, each time value is sampled, sample it from interpolated downsacled z value
    // loop the same one as tho one when templates were formed to be compliant
    // take image values first (img_vals array)

    float ag = 0;

    // the way values are stored has to be the same as with the templates tt

    if (is2D) {

        // one orthogonal is necessary
        rotation_matrix(    1,0,0, vx, vy, vz, rot); // vz should be zero
        rotation_apply(rot, 0,1,0, ux, uy, uz); // get the othogonal vector
        wx = 0;        wy = 0;        wz = 0;

        for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*1
            for (int uu = -U2; uu <= U2; ++uu) {
                //for (int ww = -W2; ww <= W2; ++ww) { // thos one will

                    int v1 = vv + V2; // 0-V
                    int u1 = uu + U2; // 0-U

                    float xcoord = vv*vx + uu*ux; // x components of all three orthogonals
                    float ycoord = vv*vy + uu*uy;

                    if ( floor(x+xcoord)<0 || ceil(x+xcoord)>img_w-1) return 0;
                    if ( floor(y+ycoord)<0 || ceil(y+ycoord)>img_h-1) return 0;

//                    printf("gginterp %f %f %f \n", x, y, z);
                    float value = interp(x+xcoord, y+ycoord, 0, img, img_w, img_h, img_l);
                    img_vals[v1*U+u1] = value;
                    ag += value;

                //}

            }
        }

    }
    else {
        // two orthogonals and differnt indexing in 2d array (the same as in template)
        rotation_matrix(    0,0,1, vx, vy, vz, rot);
        rotation_apply(rot, 1,0,0, ux, uy, uz); // get the othogonal vectors
        rotation_apply(rot, 0,1,0, wx, wy, wz);

        for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*W
            for (int uu = -U2; uu <= U2; ++uu) {
                for (int ww = -W2; ww <= W2; ++ww) { // thos one will

                    int v1 = vv + V2; // 0-V
                    int u1 = uu + U2; // 0-U
                    int w1 = ww + W2; // 0-W

                    float xcoord = vv*vx + uu*ux + ww*wx; // x components of all three orthogonals
                    float ycoord = vv*vy + uu*uy + ww*wy;
                    float zcoord = vv*vz + uu*uz + ww*wz;
                    // z coord is downscaled when sampling the image value for measurement
                    zcoord = zcoord/zDist;

                    if ( floor(x+xcoord)<0 || ceil(x+xcoord)>img_w-1) return 0;
                    if ( floor(y+ycoord)<0 || ceil(y+ycoord)>img_h-1) return 0;
                    if ( floor(z+zcoord)<0 || ceil(z+zcoord)>img_l-1) return 0;

                    float value = interp(x+xcoord, y+ycoord, z+zcoord, img, img_w, img_h, img_l);
                    img_vals[v1*U*W+w1*U+u1] = value;
                    ag += value;
                }
            }
        }

    }

    ag = ag / (U*V*W);

    // find correlation with corresponding template
    // calculate zncc, use template that corresponds to this sigma
    float corra = 0;
    float corrb = 0;
    float corrc = 0;

    for (int k = 0; k <U*W*V; k++) {
        corra += (img_vals[k]-ag) * (tt[gcsstd_idx][k]-tta[gcsstd_idx]);
        corrb += pow(img_vals[k]-ag, 2);
        corrc += pow(tt[gcsstd_idx][k]-tta[gcsstd_idx], 2);
    }

    //float corr_val = corra / sqrt(corrb*corrc);

    float corr_val = interp(x, y, z, img, img_w, img_h, img_l); // !!!
    return (corr_val<0)?0:corr_val;

}

void Tracer::save(V3DPluginCallback2 &callback, QString inimg_path) {

    unsigned char * tplt;
    tplt = new unsigned char[V*U*W];

    for (int ti = 0; ti < gcsstd_nr; ++ti) {

        QString out_str = inimg_path + QString("_template%1.tif").arg(ti);
        printf("saving %s\n",out_str.toStdString().c_str());

        V3DLONG in_sz[4];

        if (is2D) { // indexing is differnt in 2d versus 3d, v is y axis in 2d and z axis in 3d
            in_sz[0] = (long) U;
            in_sz[1] = (long) V;
            in_sz[2] = (long) W;
            in_sz[3] = (long) 1;
        }
        else {
            in_sz[0] = (long) U;
            in_sz[1] = (long) W;
            in_sz[2] = (long) V;
            in_sz[3] = (long) 1;
        }

        // save image, convert each kernel to unsigned char * with elements 0-255
        for (int kk = 0; kk < V*U*W; ++kk) {
            tplt[kk] = tt[ti][kk]*255.0;
        }

        simple_saveimage_wrapper(callback, out_str.toStdString().c_str(), tplt, in_sz, V3D_UINT8);

//        Image4DSimple outimg1;
//        outimg1.setData(knl, in_sz[0], in_sz[1], in_sz[2], in_sz[3], V3D_UINT8);
//        callback.saveImage(&outimg1, (char *)out_str.toStdString().c_str());

//        printf("what was it?\n");
//        for (int var = 0; var < Lxy*Lxy*Lz; ++var) {
//            printf("%6.2f (%d,%d,%d)", kernels[ki][var], dx[var], dy[var], dz[var]);
//        }

    }

    delete [] tplt;
    tplt = 0;

}

void Tracer::save_trace(QString inimg_path){
    // read xc, rc centroid estimations and export image_name_trace.swc
    NeuronTree curr_trace;

    for (int nodeidx = 0; nodeidx < iter_counter; ++nodeidx) {
        NeuronSWC nn;
        nn.nodeinseg_id = nodeidx+1;
        nn.type = 6;
        nn.n = nodeidx+1;
        nn.x = xc[nodeidx][0];  // gpnt_list[nodeidx][0];
        nn.y = xc[nodeidx][1];  // gpnt_list[nodeidx][1];
        nn.z = xc[nodeidx][2];  // gpnt_list[nodeidx][2];
        nn.r = rc[nodeidx];     // gpnt_list[nodeidx][3];
        nn.parent = (nodeidx==0)? -1 : nodeidx;
        curr_trace.listNeuron.append(nn);
    }

    QString swc_name = inimg_path + "_atrace.swc";
    curr_trace.name = "SingleTraceLocs";
    writeSWC_file(swc_name.toStdString().c_str(),curr_trace);
}
