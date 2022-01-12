#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "m_pd.h"
#include <mysofa.h>
#include <fftw3.h>

static t_class *mysofa_tilde_class;

#define MAX_BLOCKSIZE 16384
#define MAX_N_POINTS 3000


typedef struct _mysofa_tilde {
    t_object x_obj;
    t_float rightIR[MAX_BLOCKSIZE];
    t_float leftIR[MAX_BLOCKSIZE];
    t_float f;
    t_float spazi;
    t_float elevation;
    t_float distance;
    t_float spori;
    t_float values[4];
    t_float x,y,z,globalazi,globalori;
    t_float leftDelay;
    t_float rightDelay;
    t_float fftsize;
    t_float convsize;
    t_float nbins;
    t_float flag;
    t_symbol *filenameArg;
    t_float filter_length;
    t_float err;
    t_float delaysize;
    t_float sr;

    t_sample r_buffer[MAX_BLOCKSIZE];
    t_sample l_buffer[MAX_BLOCKSIZE];

    t_inlet *x_in2;
    t_inlet *x_in3;
    char path[2000];
    t_outlet *x_r_out;
    t_outlet *x_l_out;

    struct MYSOFA_EASY *sofa;
    struct MYSOFA_EASY *S000;
    struct MYSOFA_EASY *S015;
    struct MYSOFA_EASY *S030;
    struct MYSOFA_EASY *S045;
    struct MYSOFA_EASY *S060;
    struct MYSOFA_EASY *S075;
    struct MYSOFA_EASY *S090;
    struct MYSOFA_EASY *S105;
    struct MYSOFA_EASY *S120;
    struct MYSOFA_EASY *S135;
    struct MYSOFA_EASY *S150;
    struct MYSOFA_EASY *S165;
    struct MYSOFA_EASY *S180;
    char filename[1000];
    float *s_in; //s_in
    float *l_ir, *r_ir; //l_ir, r_ir;
    float *l_out, *r_out;//l_out, r_out;

    fftwf_complex *S_out;//S_out
    fftwf_complex *L_ir,*R_ir;//L_ir, R_ir;
    fftwf_complex *L_out, *R_out;//L_out,R_out;


    fftwf_plan plan1, plan2, plan3, plan4,plan5;

} t_mysofa_tilde;

// MAX_BLOCKSIZE 8192
//blockScale = MAX_BLOCKSIZE / blocksize;
//scaledBlocksize = blocksize * blockScale;
// blocksizeDelta = MAX_BLOCKSIZE -1 - scaledBlocksize



t_int *mysofa_tilde_perform(t_int *w) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)(w[1]);
    t_sample  *in =    (t_sample *)(w[2]);
    t_sample  *r_out =    (t_sample *)(w[3]);
    t_sample  *l_out =    (t_sample *)(w[4]);
    int       n =           (int)(w[5]);
    
    float global_Centerori,global_Centerori_by15;//Angle at which the speaker faces the listener (global coordinates)
    float localazi,localori;//Position and angle of the speaker (local coordinates)
    float globalazi_by5, globalori_by15;//Position and angle of the speaker (local coordinates)
    float localori_by15;//Listener orientation (local coordinates)
    int selectSOFA;
    float values[2];
    
    if(x->err==0){

        int i = 0;
        int j = 0;
        int l = 0;
        float realD,imagD,realS,imagS;
        float mux = 1.0/x->fftsize;
        x->nbins = x->fftsize/2 + 1;
        
        values[0] = x->spazi;
        values[1] = x->spori;
        /*
        x->values[0] = x->liori;//x->liori;
        x->values[1] = x->elevation;//x->elevation;
        x->values[2] = x->distance;//x->distance;
        x->values[3] = x->spori;//x->spori;
  */
        //get leftIR and rightIR
        if(x->globalazi != values[0] || x->globalori != values[1]){

            x->globalazi = values[0];
            //x->y = x->values[1];
            //x->z = x->values[2];
            x->globalori = values[1];
            
            globalazi_by5 = (x->globalazi)/5;
            globalazi_by5 = (int)globalazi_by5 * 5;
            if(x->globalazi > globalazi_by5 + 2.5) globalazi_by5 = globalazi_by5 + 5;
            if(globalazi_by5 == 180) localazi = 0;
            else if(globalazi_by5 < 180) localazi = 180 - globalazi_by5;
            else localazi = 540 - globalazi_by5;
            //SOFA
            if(x->globalazi < 180)global_Centerori = x->globalazi + 180;
            else global_Centerori = x->globalazi - 180;
            global_Centerori_by15 = global_Centerori/15;
            global_Centerori_by15 =  (int)global_Centerori_by15*15;
            if(global_Centerori > global_Centerori_by15 + 7.5) global_Centerori_by15 = global_Centerori_by15 + 15;
            
            globalori_by15 = (x->globalori)/15;
            globalori_by15 = (int)globalori_by15 * 15;
            if(x->globalori > globalori_by15 + 7.5) globalori_by15 = globalori_by15 + 15;
            
            localori = globalori_by15 - global_Centerori_by15;
            localori_by15 = localori/15;
            localori_by15 = (int)localori_by15*15;
            if(localori > localori_by15 + 7.5) localori_by15 = localori_by15 + 15;
            if(localori_by15 < 0) localori_by15 = localori_by15 + 360;
            //localori_by15 = 360 - localori_by15;
            if(localori_by15 == 360) localori_by15 = 0;
            //if(Sorientation > 180)Sorientation = (-1) * (Sorientation - 180);
            post("Global: Sazimuth is %f->%d, Sorientation %f->%d,SorientationToCenter is %f->%d",x->globalazi,(int)globalazi_by5,x->globalori,(int)globalori_by15,global_Centerori,(int)global_Centerori_by15);
            post("Local: Listener orientation is %d, Speaker orientation is %d",(int)localazi,(int)localori_by15);
            
            //SOFA get
            if(localori_by15 > 180) selectSOFA = 360 - (int)localori_by15;
            else selectSOFA = (int)localori_by15;
            selectSOFA = 180 - selectSOFA;
            switch(selectSOFA){
                case 0:
                    x->sofa = x->S000;
                    break;
                case 15:
                    x->sofa = x->S015;
                    break;
                case 30:
                    x->sofa = x->S030;
                    break;
                case 45:
                    x->sofa = x->S045;
                    break;
                case 60:
                    x->sofa = x->S060;
                    break;
                case 75:
                    x->sofa = x->S075;
                    break;
                case 90:
                    x->sofa = x->S090;
                    break;
                case 105:
                    x->sofa = x->S105;
                    break;
                case 120:
                    x->sofa = x->S120;
                    break;
                case 135:
                    x->sofa = x->S135;
                    break;
                case 150:
                    x->sofa = x->S150;
                    break;
                case 165:
                    x->sofa = x->S165;
                    break;
                case 180:
                    x->sofa = x->S180;
                    break;
                default:
                    error("SOFA file is nothing.");
                    break;
            }
            post("SOFA file is S%03d loaded.", selectSOFA);
            //
            x->values[0] = -localazi + 180;//x->liori;
            x->values[1] = x->elevation;//x->elevation;
            x->values[2] = x->distance;//x->distance;
            mysofa_s2c(x->values);
            x->x = x->values[0];
            x->y = x->values[1];
            x->z = x->values[2];
            post("%f,%f,%f",x->x,x->y,x->z);
            
            mysofa_getfilter_float(x->sofa,x->x,x->y,x->z,x->leftIR,x->rightIR,&x->leftDelay,&x->rightDelay);
            x->delaysize = x->rightDelay + x->leftDelay + x->fftsize;
        
            //SOFA close
            //mysofa_close_cached(x->sofa);
        }
        
        for(i = 0; i<x->fftsize ; i++){
            //signal
            if(i<n){
                x->s_in[i] = in[i];
            }
            else{
                x->s_in[i] = 0.0;
            }

            //ir
            if(i<x->filter_length){
                x->l_ir[i] = x->leftIR[i];
                x->r_ir[i] = x->rightIR[i];

            }
            else{
                x->l_ir[i] = 0.0;
                x->r_ir[i] = 0.0;
            }
        }

        //fft
        fftwf_execute(x->plan1);
        fftwf_execute(x->plan2);
        fftwf_execute(x->plan3);


        for( i = 0; i < x->nbins; i++ ){
            //left convolution
            realD = x->S_out[i][0];
            imagD = x->S_out[i][1];
            realS = x->L_ir[i][0];
            imagS = x->L_ir[i][1];
            x->L_out[i][0] = (realD * realS - imagD * imagS)*mux;
            x->L_out[i][1] = (realD * imagS + imagD * realS)*mux;

            //right convolution
            realD = x->S_out[i][0];
            imagD = x->S_out[i][1];
            realS = x->R_ir[i][0];
            imagS = x->R_ir[i][1];
            x->R_out[i][0] = (realD * realS - imagD * imagS)*mux;
            x->R_out[i][1] = (realD * imagS + imagD * realS)*mux;
        }


        //inverse fft
        fftwf_execute(x->plan4);
        fftwf_execute(x->plan5);

        //delay panning

        j=0;
        l=0;

        for(i=0; i<x->leftDelay + x->fftsize ;i++){
            if(i<x->leftDelay){
            x->l_buffer[i] = x->l_buffer[i] + 0;
            }
            else{
            x->l_buffer[i] = x->l_buffer[i] + x->l_out[i];
            }
        }
            
        
        for(i=0; i<x->rightDelay+x->fftsize;i++){
            if(i<x->rightDelay){
            x->r_buffer[i] = x->r_buffer[i] + 0;
            }
            else{
            x->r_buffer[i] = x->r_buffer[i] + x->r_out[i];
            }
        }


        //output and storage buffer
        for(i=0;i<x->fftsize;i++){
            if(i<n){
                r_out[i] = x->r_buffer[i];
                l_out[i] = x->l_buffer[i];

            }
            x->r_buffer[i] = x->r_buffer[i + n];
            x->l_buffer[i] = x->l_buffer[i + n];
         }
    }

    //r_out is right
    //l_out is left

    return (w+6);

}



void mysofa_tilde_dsp(t_mysofa_tilde *x, t_signal **sp) {
    dsp_add(mysofa_tilde_perform,
            5,
            x,
            sp[0]->s_vec, //in_signal
            sp[1]->s_vec, //out_signal_r
            sp[2]->s_vec, //out_signal_l
            sp[0]->s_n);

    int filter_length, err;
    int i=0;
    int size[8] = {128, 256, 512, 1024, 2048, 4096, 8192,16384};

    x->err = 100.0;
    x->sr = sp[0]->s_sr;

    //SOFA open
    for(int strori = 0; strori <= 180; strori = strori + 15){
            char file[2000] ="";
            char str[8] ="";
            
            strcpy(file,x->path);
            strcat(file,"/MySOFA/");
            sprintf(str, "S%03d", strori);
            strcat(file,str);
            strcat(file,"_sofa.sofa");
           
            if(strori == 0) x->S000 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 15) x->S015 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 30) x->S030 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 45) x->S045 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 60) x->S060 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 75) x->S075 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 90) x->S090 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 105) x->S105 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 120) x->S120 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 135) x->S135 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 150) x->S150 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 165) x->S165 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else if(strori == 180) x->S180 = mysofa_open_cached(file, x->sr, &filter_length, &err);
            else {
                post("S%03d SOFA file is nothing.",strori);
                break;
            }
        post("SOFA file is %s loaded.",file);

        }
    //strcat(file, "/Users/sakuraiyuki/Documents/Pd/kenkyu/sakurai-Pure-data-object-master/sakurai/Sakurai-SOFA-object-for-Pure-data-main/mysofa~/mit_kemar_normal_pinna.sofa");
    //strcat(file,"/Users/sakuraiyuki/Documents/Pd/kenkyu/sakurai-Pure-data-object-master/sakurai/Sakurai-SOFA-object-for-Pure-data-main/newMySOFA/S000_sofa.sofa");

    //x->sofa = mysofa_open(file, x->sr, &filter_length, &err);
    //mysofa_tilde_open(x, x->filenameArg);
    x->sofa = x->S000;
    x->filter_length = filter_length;
    x->convsize = x->filter_length + sp[0]->s_n - 1;
    x->err = err;
    x->spazi = 0;
    x->elevation = 0;
    x->distance = 100;
    x->spori = 0;

    if(x->err != 0){
        error( "The file could not be read.");
    }


    else{
        while(1){
            if(i==8){
                post("blocksize is too large");
                break;
            }
            if(x->convsize <= size[i]){
                x->fftsize = size[i];
                break;
            }
            i++;
        }

        post("filter_length : %f",x->filter_length);
        post("convsize : %f",x->convsize);
        post("fftsize : %f",x->fftsize);


        x->s_in = fftwf_alloc_real(x->fftsize);
        x->S_out = fftwf_alloc_complex(x->fftsize);
        x->plan1 = fftwf_plan_dft_r2c_1d(x->fftsize, x->s_in, x->S_out,FFTW_ESTIMATE);

        x->l_ir = fftwf_alloc_real(x->fftsize);
        x->L_ir = fftwf_alloc_complex(x->fftsize);
        x->plan2 = fftwf_plan_dft_r2c_1d(x->fftsize, x->l_ir, x->L_ir, FFTW_ESTIMATE);

        x->r_ir = fftwf_alloc_real(x->fftsize);
        x->R_ir = fftwf_alloc_complex(x->fftsize);
        x->plan3 = fftwf_plan_dft_r2c_1d(x->fftsize, x->r_ir, x->R_ir, FFTW_ESTIMATE);


        x->L_out = fftwf_alloc_complex(x->fftsize);
        x->l_out= fftwf_alloc_real(x->fftsize);
        x->plan4 = fftwf_plan_dft_c2r_1d(x->fftsize,x->L_out,x->l_out, FFTW_ESTIMATE);

        x->R_out = fftwf_alloc_complex(x->fftsize);
        x->r_out= fftwf_alloc_real(x->fftsize);
        x->plan5 = fftwf_plan_dft_c2r_1d(x->fftsize,x->R_out,x->r_out, FFTW_ESTIMATE);


        for(i = 0; i<x->fftsize + 10000; i++){
            x->r_buffer[i] = 0.0;
            x->l_buffer[i] = 0.0;
        }

    }

}



void mysofa_tilde_symbol(t_mysofa_tilde *x, t_symbol *s){
    x->filenameArg = s;
}





void mysofa_tilde_free(t_mysofa_tilde *x) {
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    outlet_free(x->x_r_out);
    outlet_free(x->x_l_out);

    fftwf_destroy_plan(x->plan1);
    fftwf_free(x->s_in);
    fftwf_free(x->S_out);
    fftwf_destroy_plan(x->plan2);
    fftwf_destroy_plan(x->plan3);
    fftwf_free(x->l_ir);
    fftwf_free(x->r_ir);
    fftwf_free(x->L_ir);
    fftwf_free(x->R_ir);
    fftwf_destroy_plan(x->plan4);
    fftwf_destroy_plan(x->plan5);
    fftwf_free(x->L_out);
    fftwf_free(x->R_out);
    fftwf_free(x->l_out);
    fftwf_free(x->r_out);
    mysofa_close(x->sofa);
    mysofa_cache_release_all();
}




void *mysofa_tilde_new(void) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)pd_new(mysofa_tilde_class);
    x->x_in2 = floatinlet_new(&x->x_obj, &x->spazi);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->spori);
  

    x->x_r_out = outlet_new(&x->x_obj, &s_signal);
    x->x_l_out = outlet_new(&x->x_obj, &s_signal);
    x->filenameArg =gensym(" ");
    x->err = 100.0;

    strcat(x->path,canvas_getdir(canvas_getcurrent())->s_name); /// The files should be in the same directory

    return (void *)x;
}



void mysofa_tilde_setup(void) {

    mysofa_tilde_class=class_new(gensym("mysofa~"),
                                 (t_newmethod)mysofa_tilde_new,
                                 (t_method)mysofa_tilde_free,
                                 sizeof(t_mysofa_tilde),
                                 CLASS_DEFAULT,
                                 0);

    class_addmethod(mysofa_tilde_class,
                    (t_method)mysofa_tilde_dsp,
                    gensym("dsp"),A_CANT,0);

    class_addsymbol(mysofa_tilde_class, mysofa_tilde_symbol);
    CLASS_MAINSIGNALIN(mysofa_tilde_class, t_mysofa_tilde, f);

}
