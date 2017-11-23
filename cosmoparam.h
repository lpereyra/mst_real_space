#ifndef COSMOPARAM_H
#define COSMOPARAM_H

#define RHOCRIT 2.77525E11   /* Densidad crítica del Universo [Msol h² / Mpc³] */
#define GCONS 6.67300E-20    /* Constante de Gravitación [km³ / kg / seg²]     */
#define Msol 1.9891E30       /* Masa del Sol [Kg]                              */
#define Kpc 3.08568025E16    /* Kiloparsec -> Kilometro                        */  

struct cosmoparam
{
    double  omegam			;  /* Omega Materia                         */
    double  omegal			;  /* Omega Lambda                          */
    double  omegak			;  /* Omega Curvatura                       */
		double  h0          ;  /* Parámetro de Hubble                   */
		double  lbox        ;  /* Lado del box [Kpc / h]                */
		double  vol         ;  /* Volumen del box en [Mpc**3 / h]       */
		int     npart       ;  /* Número de partículas                  */
};

extern struct cosmoparam cp;

#endif
