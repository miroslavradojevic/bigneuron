/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.

IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.

IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 * NeuronFloater_plugin.cpp
 * 2015-5-11 : by Miroslav
 */
 
#include "v3d_message.h"
#include <vector>
#include "basic_surf_objs.h"
#include "nf_dialog.h"
#include "toolbox.h"
#include "model.h"
#include "node.h"
#include "tracer.h"

#include "NeuronFloater_plugin.h"
Q_EXPORT_PLUGIN2(NeuronFloater, NeuronFloater);

using namespace std;

// default plugin parameter values
static long     channel = 1;
static int      scal    = 10;
static int      perc    = 15;
static float    znccTh  = 0.75;
static int      Ndir    = 30;
static float    angSig  = 60;
static float    gcsSig  = 3;
static int      Ni      = 100;
static int      Ns      = 50;
static float    zDist   = 2.0;
static int      saveMidres = 0;

struct input_PARA
{
    QString inimg_file;
    V3DLONG channel;

    int     scal;       // scale
    int     perc;       // foreground percentile on the scale from 1-20 (5%-100%)
    float   znccTh;     // correlation threshold
    int     Ndir;       // number of directions
    float   angSig;     // angular deviation
    float   gcsSig;     // cross-section deviation
    int     Ni;         // number of iterations
    int     Ns;         // number of states
    float   zDist;      // the distance between layers in pixels
    int     saveMidres; // save midresults

};

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA, bool bmenu);
 
QStringList NeuronFloater::menulist() const
{
	return QStringList() 
		<<tr("nf_menu")
		<<tr("about");
}

QStringList NeuronFloater::funclist() const
{
	return QStringList()
		<<tr("nf_func")
		<<tr("help");
}

const QString title = QObject::tr("NeuronFloater Plugin");

void NeuronFloater::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("nf_menu"))
	{
        bool bmenu = true;
        input_PARA PARA;

        // pick the default params
        PARA.channel = channel;
        PARA.scal = scal;
        PARA.perc = perc;
        PARA.znccTh = znccTh;
        PARA.Ndir   = Ndir;
        PARA.angSig = angSig;
        PARA.gcsSig = gcsSig;
        PARA.Ni = Ni;
        PARA.Ns = Ns;
        PARA.zDist = zDist;
        PARA.saveMidres = saveMidres;

        /*
         *  input through the menu, assign PARA fields
         */

        // legend
        vector<string> items;
        items.push_back("Channel");
        items.push_back("Scale (5+) [pix]");
        items.push_back("Foreground percentile (5-20)");
        items.push_back("Correlation (0.5-0.99)");
        items.push_back("# directions (5-100)");
        items.push_back("Angular sigma (1-360) [deg]");
        items.push_back("Gaussian Cross Section sigma (1-10) [pix]");
        items.push_back("# iterations (10-300)");
        items.push_back("# states (30-100)");
        items.push_back("Z layer dist (1-10) [pix]");
        items.push_back("Save midresults (0,1)");

        // initialization
        vector<string> inits;
        inits.push_back(QString::number(PARA.channel).toStdString().c_str());
        inits.push_back(QString::number(PARA.scal).toStdString().c_str());
        inits.push_back(QString::number(PARA.perc).toStdString().c_str());
        inits.push_back(QString::number(PARA.znccTh).toStdString().c_str());
        inits.push_back(QString::number(PARA.Ndir).toStdString().c_str());
        inits.push_back(QString::number(PARA.angSig).toStdString().c_str());
        inits.push_back(QString::number(PARA.gcsSig).toStdString().c_str());
        inits.push_back(QString::number(PARA.Ni).toStdString().c_str());
        inits.push_back(QString::number(PARA.Ns).toStdString().c_str());
        inits.push_back(QString::number(PARA.zDist).toStdString().c_str());
        inits.push_back(QString::number(PARA.saveMidres).toStdString().c_str());

//        CommonDialog dialog(items);
        CommonDialog dialog(items, inits);

        dialog.setWindowTitle(title);
        if(dialog.exec() != QDialog::Accepted) return;

        dialog.get_num("Channel", PARA.channel);
        dialog.get_num("Scale (5+) [pix]", PARA.scal);
        dialog.get_num("Foreground percentile (5-20)", PARA.perc);
        dialog.get_num("Correlation (0.5-0.99)", PARA.znccTh);
        dialog.get_num("# directions (5-100)", PARA.Ndir);
        dialog.get_num("Angular sigma (1-360) [deg]", PARA.angSig);
        dialog.get_num("Gaussian Cross Section sigma (1-10) [pix]", PARA.gcsSig);
        dialog.get_num("# iterations (10-300)", PARA.Ni);
        dialog.get_num("# states (10-100)", PARA.Ns);
        dialog.get_num("Z layer dist (1-10) [pix]", PARA.zDist);
        dialog.get_num("Save midresults (0,1)", PARA.saveMidres);

        // check input
        if(PARA.channel <= 0)                           {v3d_msg(QObject::tr("Channel is out of range")); return;}
        if(PARA.scal <= 5)                              {v3d_msg(QObject::tr("Scale is out of range")); return;;}
        if(PARA.perc < 5        || PARA.perc > 20)      {v3d_msg(QObject::tr("Percentile is out of range")); return;}
        if(PARA.znccTh < 0.5    || PARA.znccTh >= 1.0)  {v3d_msg(QObject::tr("Correlation is out of range")); return;}
        if(PARA.Ndir<5          || PARA.Ndir>100)       {v3d_msg(QObject::tr("# directions is out of range")); return;}
        if(PARA.angSig < 1      || PARA.angSig>360)     {v3d_msg(QObject::tr("Angular sigma is out of range")); return;}
        if(PARA.gcsSig < 1      || PARA.gcsSig>10)      {v3d_msg(QObject::tr("Gaussian Cross Section sigma is out of range")); return;}
        if(PARA.Ni < 10         || PARA.Ni>300      )   {v3d_msg(QObject::tr("# iterations is out of range")); return;}
        if(PARA.Ns < 10         || PARA.Ns>100      )   {v3d_msg(QObject::tr("# states is out of range")); return;}
        if(PARA.zDist < 1       || PARA.zDist>10    )   {v3d_msg(QObject::tr("Z layer dist is out of range")); return;}
        if(PARA.saveMidres<0    || PARA.saveMidres>1)   {v3d_msg(QObject::tr("saveMidres has to be 0 or 1")); return;}

        reconstruction_func(callback,parent,PARA,bmenu);

	}
	else
	{
        v3d_msg(tr("Neuron reconstruction plugin "
			"Developed by Miroslav, 2015-5-11"));
	}
}

bool NeuronFloater::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent) {
	if (func_name == tr("nf_func"))
	{
        bool bmenu = false;
        input_PARA PARA;

        vector<char*> * pinfiles = (input.size() >= 1) ? (vector<char*> *) input[0].p : 0;
        vector<char*> * pparas = (input.size() >= 2) ? (vector<char*> *) input[1].p : 0;
        vector<char*> infiles = (pinfiles != 0) ? * pinfiles : vector<char*>();
        vector<char*> paras = (pparas != 0) ? * pparas : vector<char*>();

        if(infiles.empty())     {
            fprintf (stderr, "Need input image. \n");
            return false;
        }
        else
            PARA.inimg_file = infiles[0];

        // constrain number of input parameters
        if (paras.size()!=11) {
            fprintf (stderr, "\n\nNeeds %d input parameters.\n\n", 11);
            return false;
        }

        int k=0;
        PARA.channel    = (paras.size() >= k+1)   ? atoi(paras[k])              : channel;   k++;
        PARA.scal       = (paras.size() >= k+1)   ? atoi(paras[k])              : scal;      k++;
        PARA.perc       = (paras.size() >= k+1)   ? atoi(paras[k])              : perc;      k++;
        PARA.znccTh     = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : znccTh;    k++;
        PARA.Ndir       = (paras.size() >= k+1)   ? atoi(paras[k])              : Ndir;      k++;
        PARA.angSig     = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : angSig;    k++;
        PARA.gcsSig     = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : gcsSig;    k++;
        PARA.Ni         = (paras.size() >= k+1)   ? atoi(paras[k])              : Ni;        k++;
        PARA.Ns         = (paras.size() >= k+1)   ? atoi(paras[k])              : Ns;        k++;
        PARA.zDist      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : zDist;     k++;
        PARA.saveMidres = (paras.size() >= k+1)   ? atoi(paras[k])              : saveMidres;k++;

        // check input
        if(PARA.channel <= 0)                           {v3d_msg(QObject::tr("Channel is out of range")); return false;}
        if(PARA.scal <= 5)                              {v3d_msg(QObject::tr("Scale is out of range")); return false;}
        if(PARA.perc < 5        || PARA.perc > 20)      {v3d_msg(QObject::tr("Percentile is out of range")); return false;}
        if(PARA.znccTh < 0.5    || PARA.znccTh >= 1.0)  {v3d_msg(QObject::tr("Correlation is out of range")); return false;}
        if(PARA.Ndir<5          || PARA.Ndir>100)       {v3d_msg(QObject::tr("# directions is out of range")); return false;}
        if(PARA.angSig < 1      || PARA.angSig>360)     {v3d_msg(QObject::tr("Angular sigma is out of range")); return false;}
        if(PARA.gcsSig < 1      || PARA.gcsSig>10)      {v3d_msg(QObject::tr("Gaussian Cross Section sigma is out of range")); return false;}
        if(PARA.Ni < 10         || PARA.Ni>300      )   {v3d_msg(QObject::tr("# iterations is out of range")); return false;}
        if(PARA.Ns < 10         || PARA.Ns>100      )   {v3d_msg(QObject::tr("# states is out of range")); return false;}
        if(PARA.zDist < 1       || PARA.zDist>10    )   {v3d_msg(QObject::tr("Z layer dist is out of range")); return false;}
        if(PARA.saveMidres<0    || PARA.saveMidres>1)   {v3d_msg(QObject::tr("saveMidres has to be 0 or 1")); return false;}

        reconstruction_func(callback,parent,PARA,bmenu);

	}
    else if (func_name == tr("help"))
    {

        printf("\n\n**** Usage of NeuronFloater tracing ****\n\n");
        printf("vaa3d -x NeuronFloater -f nf_func -i <inimg_file> -p <channel> <scal perc znccTh angSig gcsSig Ni Ns zDist saveMidres>\n");
        printf("inimg_file          The input image\n");
        printf("channel             Data channel for tracing. Start from 1 (default 1).\n");
        printf("scal                Scale (5+) [pix].\n"); // todo - limit the scale max value
        printf("perc                Percentile (for extracting foreground).\n"); // todo percentile 0-100 instead of 0-20
        printf("znccTh              Correlation threshold (couldn't get atof() output).\n");
        printf("Ndir                # directions (10-100)");
        printf("angSig              Angular sigma (1-360) [deg].\n");
        printf("gcsSig              Gaussian Cross Section sigma (1-10) [pix].\n");
        printf("Ni                  # iterations (10-300).\n");
        printf("Ns                  # states (10-100).\n");
        printf("zDist               Z layer dist (1-10) [pix].\n");
        printf("saveMidres          Save midresults (0,1).\n");
        printf("outswc_file         Will be named automatically based on the input image file name, so you don't have to specify it.\n\n");

	}
	else return false;

	return true;
}

bool no_overlap(int atx, int aty, int atz, int lim_xy, int lim_z, int* input_gpnt_map, int width, int height, int length) {

    bool not_ovlping = true;

    for (int xnn = atx-lim_xy; xnn <= atx+lim_xy; ++xnn) {
        for (int ynn = aty-lim_xy; ynn <= aty+lim_xy; ++ynn) {
            for (int znn = atz-lim_z; znn <= atz+lim_z; ++znn) {
                if (xnn>=0 && xnn<width  && ynn>=0 && ynn<height && znn>=0 && znn<length && input_gpnt_map[znn*width*height+ynn*width+xnn]>0) {
                    not_ovlping = false; // it overlaps with some other region
                    return not_ovlping;
                }
            }
        }
    }

    return not_ovlping;
}

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA, bool bmenu) {

    unsigned char* data1d = 0;
    V3DLONG N,M,P,sc,c;
    V3DLONG in_sz[4];
    if(bmenu)
    {
        v3dhandle curwin = callback.currentImageWindow();
        if (!curwin)
        {
            QMessageBox::information(0, "", "You don't have any image open in the main window.");
            return;
        }

        Image4DSimple* p4DImage = callback.getImage(curwin);

        if (!p4DImage)
        {
            QMessageBox::information(0, "", "The image pointer is invalid. Ensure your data is valid and try again!");
            return;
        }

        data1d = p4DImage->getRawData();
        N = p4DImage->getXDim();
        M = p4DImage->getYDim();
        P = p4DImage->getZDim();
        sc = p4DImage->getCDim();

        bool ok1;

        if(sc==1)
        {
            c=1;
            ok1=true;
        }
        else
        {
            c = QInputDialog::getInteger(parent, "Channel", "Enter channel NO:", 1, 1, sc, 1, &ok1);
        }

        if(!ok1)
            return;

        in_sz[0] = N;
        in_sz[1] = M;
        in_sz[2] = P;
        in_sz[3] = sc;

        PARA.inimg_file = p4DImage->getFileName();
    }
    else
    {
        int datatype = 0;
        if (!simple_loadimage_wrapper(callback,PARA.inimg_file.toStdString().c_str(), data1d, in_sz, datatype)) {
            fprintf (stderr, "Error happens in reading the subject file [%s]. Exit. \n",PARA.inimg_file.toStdString().c_str());
            return;
        }
        if(PARA.channel < 1 || PARA.channel > in_sz[3])
        {
            fprintf (stderr, "Invalid channel number. \n");
            return;
        }
        N = in_sz[0];
        M = in_sz[1];
        P = in_sz[2];
        sc = in_sz[3];
        c = PARA.channel;
    }

    /////////////////////////////////////////////////////////////////////////////
    v3d_msg("\n\nreconstruction_func...", 0);

    printf("initialize Model class...\n");
    Model mdl_test(PARA.scal, PARA.Ndir, P==1, PARA.zDist);


    printf("\ninitialize Tracer class...\n");
    Tracer trac_test(PARA.Ni, 50, 5, P==1, PARA.zDist, PARA.Ndir);


    // pointers to many objects (unknown) of one class
//    typedef Ant* AntPtr;
//    AntPtr * ants = new AntPtr[num_ants];
//    for (int i = 0; i < num_ants; ++i) {
//        ants[i] = new Ant();
//    }
//#include <vector>
//std::vector<Ant*> ants;
//for (int i = 0; i < num_ants; ++i) {
//    ants.push_back(new Ant());
//}
//    Node nd1(1.1, 3.4, 5.6, 7.1, Node::NOTHING); // object on stack
//    if (true) {
//        std::vector<Node> nd = std::vector<Node>();
//        nd.push_back(nd1);
//        Node * nd = new Node[3]; // object on heap
//        nd[0] = new Node(2.6, 4.5, 1.8, Node::AXON);
//        delete [] nd; nd = 0;
//    }

//    for (float xc = 0; xc <= N-1; xc+=0.42) {
//        for (float yc = 0; yc <= M-1; yc+=0.87) {
//            for (float zc = 0; zc <= P-1; zc+=2.22) {
//                float val = interp(xc, yc, zc, data1d, N, M, P);
//                if (val>0)
//                printf("[%6.1f, %6.1f, %6.1f] -> %6.1f \n", xc, yc, zc, val);
//            }
//        }
//    }

    //////////////////////////////////////////////////////////////////////////////////////////
    if (false) {
    printf("\n---descending experimental---\n");
    float * aa = new float[4]; aa[0] = 2; aa[1] = 6; aa[2] = 4; aa[3] = 8;
    int * ai = new int[4];

    printf("\nVAL (before) = \t");for (int i = 0; i < 4; ++i) {printf("%.1f\t", aa[i]);}
    printf("\nIDX (before) = \t");for (int i = 0; i < 4; ++i) {printf("%d\t", ai[i]);} // ai[i] = i;

    descending(aa, 4, ai);

    printf("\nVAL (after) = \t"); for (int i = 0; i < 4; ++i) {printf("%.1f\t", aa[i]);}
    printf("\nIDX (after) = \t"); for (int i = 0; i < 4; ++i) {printf("%d\t", ai[i]);}

    printf("\n---experimental---\n");
    }


    cout<<"----------  NeuronFloater  ----------"   <<endl;
    cout<<"channel = "  <<PARA.channel              <<endl;
    cout<<"scal = "     <<PARA.scal                 <<endl;
    cout<<"perc = "     <<PARA.perc                 <<endl;
    cout<<"znccTh = "   <<PARA.znccTh               <<endl;
    cout<<"Ndir = "     <<PARA.Ndir                 <<endl;
    cout<<"angSig = "   <<PARA.angSig               <<endl;
    cout<<"gcsSig = "   <<PARA.gcsSig               <<endl;
    cout<<"Ni = "       <<PARA.Ni                   <<endl;
    cout<<"Ns = "       <<PARA.Ns                   <<endl;
    cout<<"zDist = "    <<PARA.zDist                <<endl;
    cout<<"saveMidres = "<<PARA.saveMidres          <<endl;
    cout<<"-------------------------------------"   <<endl;

    long size = N * M * P; // N : width, M : height, P : nr. layers
    printf("\n%dx%dx%d\n", N, M, P);

    // templates (allocate before processing)
    Model mdl(PARA.scal, PARA.Ndir, P==1, PARA.zDist);
    if (PARA.saveMidres) {
        mdl.save(callback, PARA.inimg_file);
    }

    // tracer (allocate before processing)
    Tracer bayesian_tracer(PARA.Ni, PARA.Ns,    ceil(PARA.scal/2), P==1, PARA.zDist, PARA.Ndir);

    if (PARA.saveMidres) {
        bayesian_tracer.save(callback, PARA.inimg_file); // save the templates used for the measurement
        printf("\ntest trace..\n");

        // start it from the highest intensity location, just test
        int xtest = 0;
        int ytest = 0;
        int ztest = 0;
        unsigned char max_val = 0;
        for (long i = 0; i < size; ++i) {

            if (data1d[i]>max_val) {
                max_val = data1d[i];
                xtest = i % N;
                ztest = i / (N*M);
                ytest = i/N-ztest*M;
            }
        }
        // ((P==1)? 0 : (P/2.0))

//        float ang = (170/180.0)*3.14;

//        xtest = N/2;
//        ytest = M/2;
//        ztest = P/2;

        printf("\n\n--->%d|%d|%d\n", xtest, ytest, ztest);
        bayesian_tracer.trace(xtest,ytest,ztest,  1,0,0, 1.5, data1d, M, N, P, PARA.angSig, PARA.gcsSig);
        bayesian_tracer.save_trace(PARA.inimg_file); // export it to swc file
    }


    //////////////////////////////////////////////////////////////////////////////////////////
    if (true) {printf("goin' out..."); return;}
    /////////////////////////////////////////////////////////////////////////////


    float SCALE_RADIUS = 2.5; // constant that defines the scale between gaussian cross section and nerite radius

    /////////////////////////////////////////////////////////////////////////////
    unsigned char * scr = new unsigned char [size];
    for (long var = 0; var < size; ++var) scr[var] = 0;

    int range2 = scal;

    unsigned char * nbhood = new unsigned char [(2*range2+1)*(2*range2+1)];

    int cnt_nhood = 0;
    unsigned char m05 = 0;
    unsigned char m95 = 0;

    for(long i = 0; i < size; i++) {

        int x  = i % N;
        int z  = i / (N*M);
        int y  = i/N-z*M; // (M-1) - (i/N-z*M);

        if (x>=range2 && x<N-range2 && y>=range2 && y<M-range2  && data1d[i]>5 ) {

            cnt_nhood = 0;
            for (int xx = x-range2; xx <= x+range2; ++xx) {
                for (int yy = y-range2; yy <= y+range2; ++yy) {
                    nbhood[cnt_nhood] = data1d[z*N*M+yy*N+xx];
                    cnt_nhood++;
                }
            }

            m05 = quantile(nbhood, (2*range2+1)*(2*range2+1), 1,  20);
            m95 = quantile(nbhood, (2*range2+1)*(2*range2+1), 19, 20);

            scr[z*N*M+y*N+x] = m95-m05; // z*N*M+y*N+x

        }
        else
            scr[z*N*M+y*N+x] = 0; // z*N*M+y*N+x

    }

    if (PARA.saveMidres) {
        QString outimg_file = PARA.inimg_file + "_scr.tif";
        simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), scr, in_sz, V3D_UINT8);
    }

    unsigned char * scr1 = new unsigned char[size];
    for (long ii = 0; ii < size; ++ii) {scr1[ii]  = scr[ii];}
    unsigned char th = quantile(scr1, size, PARA.perc, 20);
    v3d_msg(QString("done.\tth was %1\n").arg(th), 0);
    delete [] scr1; scr1 = 0;

    unsigned char * fg = new unsigned char[size];
    for (long iii = 0; iii < size; ++iii) {fg[iii] = 0;}

    // effectively dilatation in xy
    unsigned char lmax = 0;
    bool * exam = new bool[size];                    // to bookmark the examined ones (speedup approximation)
    for (long ii = 0; ii < size; ++ii) {exam[ii] = false;}

    range2 = PARA.scal;  // largest scale chosen to dilatate

    for (long i = 0; i < size; ++i) {

        int x  = i % N;
        int z  = i / (N*M);
        int y  = i/N-z*M;

//        int iloc = z*M*N + y*N + x; // (M-1-y)
//        fg[i] = data1d[iloc];
//        printf("%d \t -> [%d, %d, %d] -> %d \t | \t %d \n", i, x, y, z, data1d[i], fg[i]);

//        if (exam[z*N*M+y*N+x]==true) { // if it was known to be close enough to the local maximum then skip it
//            continue;
//        }

        if (fg[z*N*M+y*N+x]==255) continue;

        if (x>=range2 && x<N-range2 && y>=range2 && y<M-range2) { // && scr[i]>0

//            lmax = 0;
//            for (int xx = x-range2; xx <= x+range2; ++xx) {
//                for (int yy = y-range2; yy <= y+range2; ++yy) {
//                    exam[z*N*M+yy*N+xx] = true;
//                    if(scr[z*N*M+yy*N+xx]>lmax) {
//                        lmax = scr[z*N*M+yy*N+xx];
//                    }
//                }
//            }

            if (scr[z*N*M+y*N+x]>th) { // assign to all the points
                for (int xx = x-range2; xx <= x+range2; ++xx) {
                    for (int yy = y-range2; yy <= y+range2; ++yy) {
                        if (pow(xx-x,2)+pow(yy-y,2)<=range2*range2) {
                            fg[z*N*M+yy*N+xx] = 255;
                        }
                    }
                }
            }
        }

    }

    if (PARA.saveMidres) {
        QString outimg_file = PARA.inimg_file + QString("_fg.tif");
        simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), fg, in_sz, V3D_UINT8);
    }

    // fg is calculated
    delete [] nbhood; nbhood = 0;
    delete [] scr; scr = 0;
    delete [] exam; exam = 0;

    // count foreground locs in a separate loop after all the fg[]= assignments are done
    long cnt_fg = 0;
    for (long i = 0; i < size; ++i) {if (fg[i]==255) cnt_fg++;}
    printf("\n%5.2f%\tvolume kept in fg.\n", ((float)cnt_fg/(float)(M*N*P))*100.0);

    ////////////////////////////////////////////////////////////

    // mapping i2xyz, xyz2i, i2zncc (store corr values so that they're not recalculated)
    long ** i2xyz = new long*[cnt_fg];
    for(int i = 0; i < cnt_fg; ++i)
        i2xyz[i] = new long[3];

    long * xyz2i = new long[size];
    for(long i = 0; i < size; ++i) xyz2i[i] = -1;

    cnt_fg = 0;
    for (long i = 0; i < size; ++i) {
        if (fg[i]==255) {

            int xf  = i % N;
            int zf  = i / (N*M);
            int yf  = i/N-zf*M;

            i2xyz[cnt_fg][0] = xf;
            i2xyz[cnt_fg][1] = yf;
            i2xyz[cnt_fg][2] = zf;

            xyz2i[zf*N*M+yf*N+xf] = cnt_fg;

            cnt_fg++;
        }
    }

    delete [] fg;  fg = 0;

    ////////////////////////////////////////////////////////////
    // calculate kernel corr at each fg loc.
    // if you find it (correlates above some th.), skip the surroundings of corresponding scale from calcualting
    v3d_msg("calculating corrs.", 0);

    // store calcualtions
    float * i2zncc = new float[cnt_fg];
    for (long i = 0; i < cnt_fg; ++i) i2zncc[i] = -1;

    float * i2sig = new float[cnt_fg];
    for (long i = 0; i < cnt_fg; ++i) i2sig[i] = -1;

    float * i2vx = new float[cnt_fg];
    for (long i = 0; i < cnt_fg; ++i) i2vx[i] = -1;

    float * i2vy = new float[cnt_fg];
    for (long i = 0; i < cnt_fg; ++i) i2vy[i] = -1;

    float * i2vz = new float[cnt_fg];
    for (long i = 0; i < cnt_fg; ++i) i2vz[i] = -1;

    ////////////////////////////////////////////////////////////
    std::vector<float*> gpnt_list = std::vector<float*>();  // list of gpnt locations <x,y,z,s,vx,vy,vz>
    int * gpnt_map = new int[size];                         // map that corresponds to the gpnt list
    for (long i = 0; i < size; ++i) gpnt_map[i] = -1;

    float   get_zncc, get_s, get_vx, get_vy, get_vz; // readout variables
    int     gpnt_x, gpnt_y, gpnt_z;
    float   gpnt_s, gpnt_vx, gpnt_vy, gpnt_vz;

    ////////////////////////////////////////////////////////////
    for (long i = 0; i < cnt_fg; ++i) { // go through all fg locs  cnt_fg

        int x  = i2xyz[i][0]; // get corresponding real locations in 3d
        int y  = i2xyz[i][1];
        int z  = i2xyz[i][2];

        if (gpnt_map[z*N*M+y*N+x]>=0) continue; // checked (tag==0) or assigned to some guidepoint region (tag>0)

        mdl.get_corr(x,y,z, data1d, N, M, P, get_zncc, get_s, get_vx, get_vy, get_vz);

        gpnt_map[z*N*M+y*N+x] = 0; // checked
        i2zncc[i]   = get_zncc;
        i2sig[i]    = get_s;
        i2vx[i]     = get_vx;
        i2vy[i]     = get_vy;
        i2vz[i]     = get_vz;

        if (i2zncc[i]>PARA.znccTh) {

            bool found = false;
            float lmax_zncc = i2zncc[i];

            for (int xn = x-mdl.L2xy; xn <= x+mdl.L2xy; ++xn) { // -L2xy ... +L2xy, -L2z ... +L2z neighbourhood
                for (int yn = y-mdl.L2xy; yn <= y+mdl.L2xy; ++yn) {
                    for (int zn = z-mdl.L2z; zn <= z+mdl.L2z; ++zn) {

                        if (    xn>=0 && xn<N && yn>=0 && yn<M && zn>=0 && zn<P &&
                                pow(xn-x,2)<=pow(mdl.L2xy,2) &&
                                pow(yn-y,2)<=pow(mdl.L2xy,2) &&
                                pow(zn-z,2)<=pow(mdl.L2z,2))
                        {

                            int in = xyz2i[zn*N*M+yn*N+xn];

                            if (in>=0) { // if nbr. was in the fg.

                                // get its correlation value
                                if (gpnt_map[zn*N*M+yn*N+xn] == -1) { // it was NOT calculated i2zncc[in]==-1
                                    mdl.get_corr(xn,yn,zn, data1d, N, M, P, get_zncc, get_s, get_vx, get_vy, get_vz);
                                    gpnt_map[zn*N*M+yn*N+xn] = 0;
                                    i2zncc[in]   = get_zncc;
                                    i2sig[in]    = get_s;
                                    i2vx[in]     = get_vx;
                                    i2vy[in]     = get_vy;
                                    i2vz[in]     = get_vz;
                                }
                                else {
                                    get_zncc    = i2zncc[in];
                                    get_s       = i2sig[in];
                                    get_vx      = i2vx[in];
                                    get_vy      = i2vy[in];
                                    get_vz      = i2vz[in];
                                }

                                if (get_zncc>lmax_zncc) {

                                    int A = ceil(SCALE_RADIUS*get_s);
                                    int B  = (P==1)? 0 : round(SCALE_RADIUS*(get_s/PARA.zDist));

                                    if (no_overlap(xn, yn, zn, A, B, gpnt_map, N, M, P)) {
                                        found = true;
                                        lmax_zncc   = get_zncc; // so that the next one can pick up from there

                                        gpnt_x = xn;
                                        gpnt_y = yn;
                                        gpnt_z = zn;
                                        gpnt_s    = get_s;
                                        gpnt_vx   = get_vx;
                                        gpnt_vy   = get_vy;
                                        gpnt_vz   = get_vz;
                                    }
                                }

                            }

                        }
                    } // zn
                } // yn
            } // xn

            if (found) {

                float * nn = new float[7];
                nn[0] = gpnt_x;
                nn[1] = gpnt_y;
                nn[2] = gpnt_z;
                nn[3] = gpnt_s;
                nn[4] = gpnt_vx;
                nn[5] = gpnt_vy;
                nn[6] = gpnt_vz;

                gpnt_list.push_back(nn);

                printf("%d, \t %f \t %f \t %f \t %f\n", gpnt_list.size(), gpnt_list.back()[0], gpnt_list.back()[1], gpnt_list.back()[2], gpnt_list.back()[3]);

                // label the added region in the map
                int A = ceil(SCALE_RADIUS*gpnt_s);
                int B  = (P==1)? 0 : round(SCALE_RADIUS*(gpnt_s/PARA.zDist));
                for (int xnn = gpnt_x-A; xnn <= gpnt_x+A; ++xnn) {
                    for (int ynn = gpnt_y-A; ynn <= gpnt_y+A; ++ynn) {
                        for (int znn = gpnt_z-B; znn <= gpnt_z+B; ++znn) {
                            if (xnn>=0 && xnn<N && ynn>=0 && ynn<M && znn>=0 && znn<P) {
                                gpnt_map[znn*N*M+ynn*N+xnn] = gpnt_list.size(); // >0
                            }
                        }
                    }
                }


            } // if (found)
            else { // check if the initial one (indexed with 'i') overlaps

                int A = ceil(SCALE_RADIUS*i2sig[i]);
                int B  = (P==1)? 0 : round(SCALE_RADIUS*(i2sig[i]/PARA.zDist));

                if (no_overlap(i2xyz[i][0], i2xyz[i][1], i2xyz[i][2], A, B, gpnt_map, N, M, P)) {

                    float * nn = new float[7];
                    nn[0] = i2xyz[i][0];
                    nn[1] = i2xyz[i][1];
                    nn[2] = i2xyz[i][2];
                    nn[3] = i2sig[i];
                    nn[4] = i2vx[i];
                    nn[5] = i2vy[i];
                    nn[6] = i2vz[i];

                    gpnt_list.push_back(nn); // add it to the list

                    printf("> %d, \t %f \t %f \t %f \t %f\n", gpnt_list.size(), gpnt_list.back()[0], gpnt_list.back()[1], gpnt_list.back()[2], gpnt_list.back()[3]);

                    // label the added region in the map
                    for (int xnn = i2xyz[i][0]-A; xnn <= i2xyz[i][0]+A; ++xnn) {
                        for (int ynn = i2xyz[i][1]-A; ynn <= i2xyz[i][1]+A; ++ynn) {
                            for (int znn = i2xyz[i][2]-B; znn <= i2xyz[i][2]+B; ++znn) {
                                if (xnn>=0 && xnn<N && ynn>=0 && ynn<M && znn>=0 && znn<P) {
                                    gpnt_map[znn*N*M+ynn*N+xnn] = gpnt_list.size(); // >0
                                }
                            }
                        }
                    }

                } // initial one was valid

            } // !found

        } // curr_zncc>PARA.znccTh

    } // go through all fg locs  cnt_fg

    printf("\nDOne. %d gpnt elements in list\n", gpnt_list.size());
    ////////////////////////////////////////////////////////////
    // save zncc filtering scores
    if (PARA.saveMidres) {

        unsigned char * out_zncc_plot = new unsigned char[size];

        for (long i = 0; i < cnt_fg; ++i) {

            if (i2zncc[i]>=0) {
                int x = i2xyz[i][0];
                int y = i2xyz[i][1];
                int z = i2xyz[i][2];
                out_zncc_plot[z*N*M+y*N+x] = i2zncc[i] * 255;
            }

            QString outimg_file = PARA.inimg_file + QString("_corr.tif");
            simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), out_zncc_plot, in_sz, V3D_UINT8);

            delete [] out_zncc_plot; out_zncc_plot = 0;
        }

    }

    ////////////////////////////////////////////////////////////
    delete [] i2zncc; i2zncc = 0;
    delete [] i2sig; i2sig = 0;
    delete [] i2vx; i2vx = 0;
    delete [] i2vy; i2vy = 0;
    delete [] i2vz; i2vz = 0;
    ////////////////////////////////////////////////////////////

    // save gpnts as list of isolated locations in swc
    if (PARA.saveMidres) {

        // save gpnts in form of swc file
        NeuronTree nt_locs;

        for (int nodeidx = 0; nodeidx < gpnt_list.size(); ++nodeidx) {

            NeuronSWC nn;
            nn.nodeinseg_id = nodeidx+1;
            nn.type = 6;
            nn.n = nodeidx+1;
            nn.x = gpnt_list[nodeidx][0];
            nn.y = gpnt_list[nodeidx][1];
            nn.z = gpnt_list[nodeidx][2];
            nn.r = gpnt_list[nodeidx][3];
            nn.parent = -1;
            nt_locs.listNeuron.append(nn);

        }

        QString swc_name = PARA.inimg_file + "_gpnt.swc";
        nt_locs.name = "GuidepointLocs";
        writeSWC_file(swc_name.toStdString().c_str(),nt_locs);

    }

    ////////////////////////////////////////////////////////////
    // cahsing (using the tracer instance)


    ////////////////////////////////////////////////////////////

    delete [] gpnt_map; gpnt_map = 0;
    for (int i = 0; i < cnt_fg; ++i) delete [] i2xyz[i];
    delete [] i2xyz; i2xyz = 0;
    delete [] i2sig; i2sig = 0;
    delete [] i2vx; i2vx = 0;
    delete [] i2vy; i2vy = 0;
    delete [] i2vz; i2vz = 0;
    /*
     *
     */

    //Output
//    NeuronTree nt;
//	QString swc_name = PARA.inimg_file + "_NeuronFloater.swc";
//	nt.name = "NeuronFloater";
//    writeSWC_file(swc_name.toStdString().c_str(),nt);

    if(!bmenu) { if(data1d) {delete []data1d; data1d = 0;}}

//    v3d_msg(QString("Now you can drag and drop the generated swc fle [%1] into Vaa3D.").arg(swc_name.toStdString().c_str()),bmenu);

    return;
}

