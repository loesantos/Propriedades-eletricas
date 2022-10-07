//======================= Gerador de geometrias ==============================//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#define solido  0
#define fluido  1

using namespace std;

void gravacao ( int* );
void gera_cilindro ( int*, int );

int nx;
int ny;
int nz;

//----------------------------------------------------------------------------//

int main () {
srand(23);

int *ini_meio;
int *meio;

//----------------------------------------------------------------------------//

cout << "\n\nDimention X:    ";

cin >> nx;

//----------------------------------------------------------------------------//

cout << "\n\nDimention Y:    ";

cin >> ny;

//----------------------------------------------------------------------------//

cout << "\n\nDimention z:    ";

cin >> nz;

//----------------------------------------------------------------------------//

int ptos_meio = nx * ny * nz;

ini_meio = (int*)calloc((ptos_meio),sizeof(int));

//----------------------------------------------------------------------------//

for (int z=0; z<nz; z++ )
    {
    for (int y=0; y<ny; y++ )
        {
        for (int x=0; x<nx; x++ )
            {
            meio = ini_meio + x + y*nx + z*nx*ny;
            *meio = fluido;
            
            //if ( y == 0 || y == ny - 1 ) *meio = solido;            
            }
        }
    }

//----------------------------------------------------------------------------//
/*/
int raio;

if ( nz < ny ) raio = ( nz - 1 ) / 2;
else raio = ( ny - 1 ) / 2;

gera_cilindro ( ini_meio, raio/ );
/*/
//----------------------------------------------------------------------------//

gravacao ( ini_meio );
}


//===================== Gravação dos Resultados ==============================//

void gravacao ( int *ini_meio )
{

int *meio;
char nome[50] ="meio.txt";

ofstream f_out(nome);

cout << "\n\nGravando." << endl << endl;

/*/
f_out << 1.0 << endl;

f_out << nx << " " << ny << " " << nz << endl;

for (int z=0; z<nz; z++ )
    {
    for (int y=0; y<ny; y++ )
        {
        for (int x=0; x<nx; x++ )
            {
            meio = ini_meio + x + y*nx + z*nx*ny;

            f_out << *meio << "     ";
            }
        f_out << endl;
        }
    }

f_out.close();
/*/

ofstream fmeio( "meio.vtk" );

fmeio << "# vtk DataFile Version 2.0" << endl;
fmeio << "Geometria" << endl;
fmeio << "ASCII" << endl;
fmeio << "DATASET STRUCTURED_POINTS" << endl;
fmeio << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
fmeio << "ASPECT_RATIO 1 1 1" << endl;
fmeio << "ORIGIN 0 0 0" << endl;
fmeio << "POINT_DATA " << nx * ny * nz << endl;
fmeio << "SCALARS Geometria float" << endl;
fmeio << "LOOKUP_TABLE default" << endl;

for ( int z = 0; z < nz; z++ )
    {
    for ( int y = 0; y < ny; y++ )
        {
        for ( int x = 0; x < nx; x++ )
            {
            meio = ini_meio + x + y*nx + z*nx*ny;

            fmeio << *meio << "     ";
            }
        fmeio << endl;
        }
    }

fmeio.close();

}

//======================== Gera cilindro =====================================//

void gera_cilindro ( int *ini_meio, int raio )
{
int *meio;
int in_cil;

int y0 = ny / 2;
int z0 = nz / 2;

for (int y = 0; y < ny; y++ )
    {
    for ( int z = 0; z < nz ;z++ )
        {
        for ( int x = 0; x < nx; x++ )
            {
            meio = ini_meio + x + y*nx + z*nx*ny;

            in_cil = (y-y0)*(y-y0) + (z-z0)*(z-z0);

             if ( in_cil < raio*raio ) *meio = fluido;
             }
         }
    }
}

//============================================================================//

