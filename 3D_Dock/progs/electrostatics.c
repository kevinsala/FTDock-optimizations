/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"
#include "xmmintrin.h"
#include "immintrin.h"

void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int	residue , atom ;

/************/

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      This_Structure.Residue[residue].Atom[atom].charge = 0.0 ;

      /* peptide backbone */

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " N  " ) == 0 ) {
        if( strcmp( This_Structure.Residue[residue].res_name , "PRO" ) == 0 ) {
          This_Structure.Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          This_Structure.Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 ) This_Structure.Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " O  " ) == 0 ) {
        This_Structure.Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == This_Structure.length  ) This_Structure.Residue[residue].Atom[atom].charge = -1.00 ;
      }

      /* charged residues */

      if( ( strcmp( This_Structure.Residue[residue].res_name , "ARG" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "ASP" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "GLU" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "LYS" ) == 0 ) && ( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NZ " ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  1.00 ;

    }
  }

/************/

}

#define ELECTRIC_FIELD(_position, _phy) 										\
{															\
  float _distance = PYTHAGORAS(coord1[_position], coord2[_position], coord3[_position], x_centre, y_centre, z_centre);	\
  _distance = (_distance < 2.0) ? 2.0 : _distance;									\
  float _epsilon = 80;													\
  _epsilon = (_distance >= 2.0 & _distance <= 6.0) ? 4 : _epsilon;							\
  _epsilon = (_distance > 6.0 & _distance < 8.0) ? (38 * _distance) - 224 : _epsilon;					\
															\
  _phy += (charge[_position] / (_epsilon * _distance));									\
}


// #define ELECTRIC_FIELD_V(_pcoord1, _pcoord2, _pcoord3, _pcharge, _phy) 					\
// {													\
//   float _distance = PYTHAGORAS(*(_pcoord1), *(_pcoord2), *(_pcoord3), x_centre, y_centre, z_centre);	\
//   _distance = (_distance < 2.0) ? 2.0 : _distance;							\
//   float _epsilon = 80;											\
//   _epsilon = (_distance >= 2.0 & _distance <= 6.0) ? 4 : _epsilon;					\
//   _epsilon = (_distance > 6.0 & _distance < 8.0) ? (38 * _distance) - 224 : _epsilon;			\
// 													\
//   _phy += (*(_pcharge) / (_epsilon * _distance));							\
// }


float _mm256_reduce_add_ps(__m256 a) {
    __m128 hi = _mm256_extractf128_ps(a, 1);
    __m128 lo = _mm256_extractf128_ps(a, 0);
    lo = _mm_add_ps(hi, lo);
    float ret;
    ret  = _mm_extract_ps(lo, 0);
    ret += _mm_extract_ps(lo, 1);
    ret += _mm_extract_ps(lo, 2);
    ret += _mm_extract_ps(lo, 3);
    return ret;
}


#define ELECTRIC_FIELD_V(_pcoord1, _pcoord2, _pcoord3, _pcharge, _phy) 		\
{										\
  __m256 _x = _mm256_load_ps(_pcoord1);						\
  __m256 _y = _mm256_load_ps(_pcoord2);						\
  __m256 _z = _mm256_load_ps(_pcoord3);						\
  __m256 _c = _mm256_load_ps(_pcharge);						\
										\
  /* x = x1 - x2; y = y1 - y2; z = z1 - z2; */					\
  _x = _mm256_sub_ps(_x, maskxcentre);						\
  _y = _mm256_sub_ps(_y, maskycentre);						\
  _z = _mm256_sub_ps(_z, maskzcentre);						\
										\
  /* x = x * x; y = y * y; z = z * z; */					\
  _x = _mm256_mul_ps(_x, _x);							\
  _y = _mm256_mul_ps(_y, _y);							\
  _z = _mm256_mul_ps(_z, _z);							\
										\
  /* xyz = x + y + z */								\
  __m256 _xyz = _mm256_add_ps(_x, _y);						\
  _xyz = _mm256_add_ps(_xyz, _z);						\
										\
  /* pyth = SQRT(xyz) */							\
  __m256 _pyth = _mm256_sqrt_ps(_xyz);						\
										\
  /* distance = MAX(pyth, 2) */							\
  __m256 _distance = _mm256_max_ps(_pyth, mask2);				\
										\
  /* _ge2_and_le6 = (distance => 2.0 & distance <= 6.0) */			\
  __m256 _ge2 = _mm256_cmp_ps(_distance, mask2, 29);				\
  __m256 _le6 = _mm256_cmp_ps(_distance, mask6, 18);				\
  __m256 _ge2_and_le6 = _mm256_and_ps(_ge2, _le6);				\
										\
  /* ge2_and_le6 = (distance > 6.0 & distance < 8.0) */				\
  __m256 _g6 = _mm256_cmp_ps(_distance, mask6, 30);				\
  __m256 _l8 = _mm256_cmp_ps(_distance, mask8, 17);				\
  __m256 _g6_and_l8 = _mm256_and_ps(_g6, _l8);					\
										\
  /* epsilon = (distance >= 2.0 & distance <= 6.0) ? 4 : 80.0 */		\
  __m256 _epsilon = _mm256_blendv_ps(mask80, mask4, _ge2_and_le6);		\
										\
  /* epsilon = (distance > 6.0 & distance < 8.0) ? (38 * distance) - 		\
   * 224 : epsilon */								\
  __m256 _38_mul_dist = _mm256_mul_ps(mask38, _distance);			\
  __m256 _38_mul_dist_sub_224 = _mm256_sub_ps(_38_mul_dist, mask224);		\
  _epsilon = _mm256_blendv_ps(_epsilon, _38_mul_dist_sub_224, _g6_and_l8);	\
										\
  /* phi += charge / (_epsilon * _distance) */					\
  __m256 _result = _mm256_mul_ps(_epsilon, _distance);				\
  _result = _mm256_div_ps(_c, _result);						\
  _phy = _mm256_add_ps(_phy, _result);						\
}


/************************/

void electric_field( struct Structure This_Structure , float grid_span , int grid_size , float *grid ) {

/************/

  /* Counters */

  int residue , atom;

  /* Co-ordinates */

  int x , y , z, atom_number;
  float x_centre , y_centre , z_centre;

  /* Variables */

  float phi;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {
        grid[gaddress(x,y,z,grid_size)] = (float)0 ;
      }
    }
  }

/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;
  
  /* Calculates the number of atoms */
  atom_number = 0;
  for (residue = 1; residue <= This_Structure.length; residue++) {
    atom_number += This_Structure.Residue[residue].size;
  }
  
  /* Allocates memory for all the needed attibutes of the atoms */
  float *coord1, *coord2, *coord3, *charge;
  int error;
  error  = posix_memalign((void **) &coord1, 16, atom_number * sizeof(float));
  error |= posix_memalign((void **) &coord2, 16, atom_number * sizeof(float));
  error |= posix_memalign((void **) &coord3, 16, atom_number * sizeof(float));
  error |= posix_memalign((void **) &charge, 16, atom_number * sizeof(float));
  
  if (error) {GENERAL_MEMORY_PROBLEM};
  
  /* Initializes the consecutive arrays */
  int atom_position = 0;
  for (residue = 1; residue <= This_Structure.length; residue++) {
    Amino_Acid * aminoacid = &(This_Structure.Residue[residue]);
    int aminoacid_size = aminoacid->size;
    
    for (atom = 1; atom <= aminoacid_size; atom++) {
      Atom * atom_obj = &(aminoacid->Atom[atom]);
      if (atom_obj->charge != 0) {
	coord1[atom_position] = atom_obj->coord[1];
	coord2[atom_position] = atom_obj->coord[2];
	coord3[atom_position] = atom_obj->coord[3];
	charge[atom_position] = atom_obj->charge;
	++atom_position;
      }
    }
  }
  
  /* Sets the actual number of computable atoms */
  atom_number = atom_position;
  
  /* Initializes constant values */
  __m256 mask2   = _mm256_set1_ps(2.0);
  __m256 mask4   = _mm256_set1_ps(4.0);
  __m256 mask6   = _mm256_set1_ps(6.0);
  __m256 mask8   = _mm256_set1_ps(8.0);
  __m256 mask80  = _mm256_set1_ps(80.0);
  __m256 mask38  = _mm256_set1_ps(38.0);
  __m256 mask224 = _mm256_set1_ps(224.0);

  for (x = 0; x < grid_size; x++) {

    printf( "." ) ;

    x_centre = gcentre(x, grid_span, grid_size);
    __m256 maskxcentre = _mm256_set1_ps(x_centre);

    for (y = 0; y < grid_size; y++) {

      y_centre = gcentre(y, grid_span, grid_size);
      __m256 maskycentre = _mm256_set1_ps(y_centre);

      for (z = 0; z < grid_size; z++) {

        z_centre = gcentre(z, grid_span, grid_size);
	__m256 maskzcentre = _mm256_set1_ps(z_centre);
	
	phi = 0.0;
	__m256 vphi = _mm256_set1_ps(0.0);
	
	float * pcoord1 = coord1;
	float * pcoord2 = coord2;
	float * pcoord3 = coord3;
	float * pcharge = charge;
	
	for (atom_position = 0; atom_position < atom_number - 7; atom_position += 8) {
	  ELECTRIC_FIELD_V(pcoord1, pcoord2, pcoord3, pcharge, vphi);
	  pcoord1+=8; pcoord2+=8; pcoord3+=8; pcharge+=8;
	}
	
	for ( ; atom_position < atom_number; ++atom_position)
	  ELECTRIC_FIELD(atom_position, phi);

        grid[gaddress(x, y, z, grid_size)] = phi + _mm256_reduce_add_ps(vphi);
      }
    }
  }

  printf( "\n" ) ;

/************/

  return ;

}



/************************/



void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , float *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	x_low , x_high , y_low , y_high , z_low , z_high ;

  float		a , b , c ;
  float		x_corner , y_corner , z_corner ;
  float		w ;

  /* Variables */

  float		one_span ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        grid[gaddress(x,y,z,grid_size)] = (float)0 ;

      }
    }
  }

/************/

  one_span = grid_span / (float)grid_size ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( This_Structure.Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( This_Structure.Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( This_Structure.Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        x_high = x_low + 1 ;
        y_high = y_low + 1 ;
        z_high = z_low + 1 ;

        a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        for( x = x_low ; x <= x_high  ; x ++ ) {
 
          x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;

          for( y = y_low ; y <= y_high  ; y ++ ) {

            y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;

            for( z = z_low ; z <= z_high  ; z ++ ) {

              z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;

              w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;

              grid[gaddress(x,y,z,grid_size)] += (float)( w * This_Structure.Residue[residue].Atom[atom].charge ) ;

            }
          }
        }

      }

    }
  }

/************/

  return ;

}



/************************/



void electric_field_zero_core( int grid_size , float *elec_grid , float *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x , y , z ;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        if( surface_grid[gaddress(x,y,z,grid_size)] == (float)internal_value ) elec_grid[gaddress(x,y,z,grid_size)] = (float)0 ;

      }
    }
  }

/************/

  return ;

}
