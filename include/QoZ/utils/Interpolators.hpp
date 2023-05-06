//
// Created by Kai Zhao on 9/1/20.
//

#ifndef SZ_INTERPOLATORS_HPP
#define SZ_INTERPOLATORS_HPP
namespace QoZ {
    template<class T>
    inline T interp_linear(T a, T b) {
        return (a + b) / 2;
    }

    template<class T>
    inline T interp_linear1(T a, T b) {
        return -0.5 * a + 1.5 * b;
    }

    template<class T>
    inline T interp_quad_1(T a, T b, T c) {
        return (3 * a + 6 * b - c) / 8;
    }

    template<class T>
    inline T interp_quad_2(T a, T b, T c) {
        return (-a + 6 * b + 3 * c) / 8;
    }

    template<class T>
    inline T interp_quad_3(T a, T b, T c) {
        return (3 * a - 10 * b + 15 * c) / 8;
    }

    template<class T>
    inline T interp_cubic(T a, T b, T c, T d) {
        return (-a + 9 * b + 9 * c - d) / 16;
        //return -0.06368435202786181*a+0.5731591682507563*b+0.5731591682507563*c-0.06368435202786181*d;
        return (-a+3*b+3*c-d)/4;
    }

    template<class T>
    inline T interp_2d(T a, T b, T c, T d) {
        return ( a + b + c + d) / 4;
    }

    template<class T>
    inline T interp_3d(T a, T b, T c, T d, T e,T f) {
        return ( a + b + c + d + e + f ) / 6;
    }

    template<class T>
    inline T lorenzo_2d(T a, T b, T c) {
        return (b+c-a);
    }
    
    template<class T>
    inline T lorenzo_3d(T a, T b, T c, T d, T e,T f,T g) {
        return (a-b-c+d-e+f+g);
    }

    template<class T>
    inline T interp_ave3(T a, T b, T c) {
        return (a + b+c) / 3;
    }

    template<class T>
    inline T interp_cubic_front(T a, T b, T c, T d) {
        return (5 * a + 15 * b - 5 * c + d) / 16;
    }

    template<class T>
    inline T interp_cubic_front_2(T a, T b, T c, T d) {
        return ( a + 6 * b - 4 * c + d) / 4;
    }

    template<class T>
    inline T interp_cubic_back_1(T a, T b, T c, T d) {
        return (a - 5 * b + 15 * c + 5 * d) / 16;
    }

    template<class T>
    inline T interp_cubic_back_2(T a, T b, T c, T d) {
        return (-5 * a + 21 * b - 35 * c + 35 * d) / 16;
    }

    template<class T>
    inline T interp_cubic2(T a, T b, T c, T d) {
        return (-3 * a + 23 * b + 23 * c - 3 * d) / 40;
    }

    template<class T>
    inline T interp_akima(T a, T b, T c, T d) {//what is it?
        T t0 = 2 * b - a - c;
        T t1 = 2 * c - b - d;
        T abt0 = fabs(t0);
        T abt1 = fabs(t1);
        if (fabs(abt0 + abt1) > 1e-9) {
            return (b + c) / 2 + (t0 * abt1 + t1 * abt0) / 8 / (abt0 + abt1);
        } else {
            return (b + c) / 2;
        }
    }

    template<class T>
    inline T interp_pchip(T a, T b, T c, T d) {//what is it?
        T pchip = (b + c) / 2;
        if ((b - a < 0) == (c - b < 0) && fabs(c - a) > 1e-9) {
            pchip += 1 / 4 * (b - a) * (c - b) / (c - a);
        }
        if ((c - b < 0) == (d - c < 0) && fabs(d - b) > 1e-9) {
            pchip -= 1 / 4 * (c - b) * (d - c) / (d - b);
        }
        return pchip;
    }

    /*
    template<class T>
    inline T lanczos(T x, int a) {
        if(x==0)
            return 0;
        else if (fabs(x)>a)
            return 1;
        else{
            T pix=M_PI*x;
            return a*sin(pix)*sin(pix/a)/(pix*pix);
        }
    }
    template<class T>
    inline T interp_lanczos_2(T a, T b, T c, T d) {

        return a*lanczos(1.5,2)+b*lanczos(0.5,2)+c*lanczos(-0.5,2)+d*lanczos(-1.5,2);
        
        -0.06368435202786181
        0.5731591682507563

        
    }
    */

}
#endif //SZ_INTERPOLATORS_HPP
