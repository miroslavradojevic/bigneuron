#ifndef TRACER_H
#define TRACER_H

#include <v3d_interface.h>

class Tracer
{
public:

    int Niterations;
    int Nstates;
    int Ndirs;
    bool is2D;
    float zDist;

    int     U, W, V, U2, W2, V2; // V = 2*V2+1 ~ [vx, vy, vz]

    float * gcsstd;
    int     gcsstd_nr;

    // template size U*W*V, W is added for 3d (it is 1 in 2d)
    // z components of each are scaled down

    float *     img_vals;   // U*W*V

    float       start_px; // each trace starts from oriented 3d point
    float       start_py;
    float       start_pz;
    float       start_vx;
    float       start_vy;
    float       start_vz;

    float **    tt;         // templates, U*W*V
    float *     tta;        // template averages

    float ***  xt;             // bayesian filtering states     Niterations*Nstates*(x,y,z,vx,vy,vz,gcsstd)
    float **   wt;             // bayesian filtering weights    Niterations*Nstates
    int **     prevt;          // state index of the precessor (val. 0-Nstates) Niterations*Nstates

    int        iter_counter;   // number of iterations carried out before stop

    float **   xc;             // centroid locations  Niterations*3
    float *    rc;             // centroid radiuses   Niterations
    int        last_queue_element_checked; // remembers at which iteration the last

    // sequential bayesian filtering
    // first prediction (p0, v0, pred_w0) will make Nstates predictions
    // all follow-up (p1, v1, pred_w1) will make Ndirs predictions
    float **    p0;      // Nstates x 3
    float **    v0;      // Nstates x 3
    float *     pred_w0; // Nstates*gcsstd_nr

    float **    p1;         // Ndirs x 3
    float **    v1;         // Ndirs x 3
    float *     pred_w1;    // Ndirs*gcsstd_nr

    // auxiliary variables assigned with values at each recursion iteration
    float **    trans_xyz0;     // 1*Nstates*gcsstd_nr x 7
    float *     pties0;         // 1*Nstates*gcsstd_nr
    float *     lhoods0;        // 1*Nstates*gcsstd_nr
    int *       sort_idx0;      // 1*Nstates*gcsstd_nr      auxilliary variable
    bool*       check0;

    float **    trans_xyz1;     // Nstates*Ndirs*gcsstd_nr x 7
    float *     pties1;         // Nstates*Ndirs*gcsstd_nr
    float *     lhoods1;        // Nstates*Ndirs*gcsstd_nr
    int *       sort_idx1;      // Nstates*Ndirs*gcsstd_nr      auxilliary variable
    bool*       check1;         // mark those that are checked (will be used for tagging the picked locations)

    float **    rot;// = new float[3][3]; // for rotations
    float       ux, uy, uz;
    float       wx, wy, wz;

    // each of the Ndirs will correspond to directional sampling in 0th direction (similar as in Model class)
    // the values will be rotated according to current direction
    float *         vx_0;                 // Nstates of such at initialization
    float *         vy_0;                 //
    float *         vz_0;                 //

    float *         vx_1;     // regular, Ndirs of such
    float *         vy_1;     //
    float *         vz_1;     //

    static float gcsstd_step;
    static float gcsstd_min;
    static float step;

    Tracer(int _Niterations, int _Nstates, int _neuron_radius, bool _is2D, float _zDist, int _Ndirs);
    ~Tracer();

    void calculate_dirs(float * _vx, float * _vy, float * _vz, bool is2D, int nr_dirs);
    void rotation_matrix(float a1, float a2, float a3, float b1, float b2, float b3, float** R);
    void rotation_apply(float** R, float v1, float v2, float v3, float &out1, float &out2, float &out3);

    float interp(float atX, float atY, float atZ, unsigned char * img, int width, int height, int length);

    void trace(float x,  float y,  float z,
               float vx, float vy, float vz,
               float gcsstd,
               unsigned char * img,
               int          img_width,
               int          img_height,
               int          img_length,
               float angula_diff_std_deg,
               float gcsstd_diff_std_pix); // trace will be stored in xc, rc

    void predict0(  float px, float py, float pz,
                    float vx, float vy, float vz,
                    float gcsstd_pv,
                    float angula_diff_std_deg, float gcsstd_diff_std_pix);

    void predict1(      float px, float py, float pz,
                        float vx, float vy, float vz,
                        float vx_prev, float vy_prev, float vz_prev,
                        float gcsstd_pv,
                        float angula_diff_std_deg, float gcsstd_diff_std_pix);

    void iter_beg(float                   xloc,
              float                   yloc,
              float                   zloc,
              float                   vxloc,
              float                   vyloc,
              float                   vzloc,
              float                   gcsstdcloc,
              unsigned char *         img,
              int                     img_width,
              int                     img_height,
              int                     img_length,
              float                   angula_diff_std_deg,
              float                   gcsstd_diff_std_pix);

    void iter(int                     iter,
              unsigned char *         img,
              int                     img_width,
              int                     img_height,
              int                     img_length,
              float                   angula_diff_std_deg,
              float                   gcsstd_diff_std_pix);

    float zncc(float x, float y, float z, float vx, float vy, float vz, int gcsstd_idx, unsigned char * img, int img_w, int img_h, int img_l);

    void save(V3DPluginCallback2 &callback, QString inimg_path);

    void save_trace(QString inimg_path);

};

#endif // TRACER_H
