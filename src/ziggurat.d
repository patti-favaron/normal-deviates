/*

Dlang implementation of ziggurat random number generator.

This is a port of John Burkardt's "ziggurat.c" and "ziggurat.h".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

Credit:

    John Burkardt,  author of a multi-language library
                    providing the same functionalities.

                    You may find his many works (including
                    ziggurat.h and ziggurat.c) here:

                        https://people.sc.fsu.edu/~jburkardt/

Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

Copyright 2023 Patrizia Favaron

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

module ziggurat;

import std.math;
import std.conv;

struct RandomDeviates {

    // Internal data for exponential generator
    uint[256]  ke;
    float[256] fe;
    float[256] we;

    // Internal data for normal generator
    uint[128]  kn;
    float[128] fn;
    float[128] wn;

    // Constructor
    this() {

        // Initialize data for exponential generator
        int i;
        double de = 7.697117470131487;
        double m2 = 2147483648.0;
        double q;
        double te = 7.697117470131487;
        double ve = 3.949659822581572E-03;

        q = ve / exp(-de);

        this.ke[0] = to!uint(( de / q ) * m2);
        this.ke[1] = 0;

        this.we[0]   = to!float(q / m2);
        this.we[255] = to!float(de / m2);

        this.fe[0]   = 1.0;
        this.fe[255] = to!float(exp(-de));

        for (i = 254; 1 <= i; i--)
        {
            de = - log(ve / de + exp(-de));
            this.ke[i+1] = to!uint((de / te) * m2);
            te = de;
            this.fe[i] = to!float(exp(-de));
            this.we[i] = to!float(de / m2);
        }

        // Initialize data for normal generator
        double dn = 3.442619855899;
        double m1 = 2147483648.0;
        double tn = 3.442619855899;
        double vn = 9.91256303526217E-03;

        q = vn / exp(-0.5 * dn * dn);

        this.kn[0] = to!uint((dn / q) * m1);
        this.kn[1] = 0;

        this.wn[0]   = to!float(q / m1);
        this.wn[127] = to!float(dn / m1);

        this.fn[0]   = 1.0;
        this.fn[127] = to!float(exp(-0.5 * dn * dn));

        for(i = 126; 1 <= i; i--)
        {
            dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn )));
            this.kn[i+1] = to!uint((dn / tn) * m1);
            tn = dn;
            this.fn[i] = to!float(exp(-0.5 * dn * dn));
            this.wn[i] = to!float(dn / m1);
        }

    }

    // Simple linear congruential generator
    uint congruential(ref uint i) {

        i = 69069 * i + 1234567;

        uint value = i;
        return value;

    }


    // KISS generator
    uint kiss_seeded(ref uint jcong, ref uint jsr, ref uint w, ref uint z) {

        uint value = pow(this.multiply_with_carry(w, z), this.congruential(jcong)) + this.shr3(jsr);
        return value;

    }


    // Multiply-with-carry generator
    uint multiply_with_carry(ref uint w, ref uint z) {

        z = 36969 * (z & 65535) + (z >> 16);
        w = 18000 * (w & 65535) + (w >> 16);

        uint value = (z << 16) + w;
        return value;

    }


    // SHR3 generator
    uint shr3(ref uint jsr) {

        uint value = jsr;

        jsr = pow(jsr, (jsr << 13));
        jsr = pow(jsr, (jsr >> 17));
        jsr = pow(jsr, (jsr <<  5));

        value += jsr;
        return value;

    }


    // Esponential generator
    float exponential(ref uint jsr) {

        uint  iz;
        uint  jz;
        float value;
        float x;

        jz = shr3(jsr);
        iz = (jz & 255);

        if (jz < this.ke[iz])
        {
            value = to!float(jz) * this.we[iz];
        }
        else
        {
            for ( ; ; )
            {
                if (iz == 0)
                {
                    value = 7.69711 - log(this.uniform(jsr));
                    break;
                }

                x = to!float(jz) * this.we[iz];

                if(this.fe[iz] + this.uniform(jsr) * (this.fe[iz-1] - this.fe[iz]) < exp(-x))
                {
                    value = x;
                    break;
                }

                jz = this.shr3(jsr);
                iz = (jz & 255);

                if(jz < this.ke[iz])
                {
                    value = to!float(jz) * this.we[iz];
                    break;
                }

            }
        }

        return value;

    }


    // Normal generator
    float normal(ref uint jsr)
    {
        int hz;
        uint iz;
        float r = 3.442620;
        float value;
        float x;
        float y;

        hz = to!int(this.shr3(jsr));
        iz = (hz & 127);

        if(fabs(hz) < this.kn[iz])
        {
            value = to!float(hz) * this.wn[iz];
        }
        else
        {
            for ( ; ; )
            {
                if (iz == 0)
                {
                    for ( ; ; )
                    {
                        x = -0.2904764 * log(this.uniform(jsr));
                        y = -log(this.uniform(jsr));
                        if (x * x <= y + y)
                        {
                            break;
                        }
                    }

                    if ( hz <= 0 )
                    {
                        value = -r - x;
                    }
                    else
                    {
                        value = +r + x;
                    }
                    break;
                }

                x = to!float(hz) * this.wn[iz];

                if(this.fn[iz] + this.uniform(jsr) * (this.fn[iz-1] - this.fn[iz]) 
                    < exp(-0.5 * x * x))
                {
                    value = x;
                    break;
                }

                hz = to!int(this.shr3(jsr));
                iz = ( hz & 127 );

                if(fabs(hz) < this.kn[iz])
                {
                    value = to!float(hz) * this.wn[iz];
                    break;
                }
            }
        }

        return value;
    }


    // Uniform generator
    float uniform(ref uint jsr)
    {
        uint jsr_input;
        float value;

        jsr_input = jsr;

        jsr = pow(jsr, (jsr << 13));
        jsr = pow(jsr, (jsr >> 17));
        jsr = pow(jsr, (jsr <<  5));

        value = fmod(0.5 + to!float(jsr_input + jsr) / 65536.0 / 65536.0, 1.0);

        return value;
    }

}
