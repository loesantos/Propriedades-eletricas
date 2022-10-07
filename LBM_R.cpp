//====================== LB --> Lattice Boltzmann ================================================//

#include <omp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#define nvel 15

#define dim  3

double W[nvel];

const double one_over_c_s2 = 3.0;

//===================== Declaração de Funções ====================================================//

//	Lê arquivos de inicialização

void read_data ( string&, double&, int&, int&, double&, double&, double& );

// Lê arquivo de geometria

int read_geo ( string, int*, int, int );

// Calcula o produto interno

double prod_int ( double*, double* );

// Arredonda um número

int sum_round ( double  );

// Define os vetores da rede D3Q15

void def_lattice_d3q15 ( double*, double* );

// Calcula a distribuição de equilíbrio

void dist_eq ( double, double, double, double, double*, double*, double*, double );

// Calcula densidade e velocidades

void calcula ( double*, double*, double&, double&, double&, double& );

// Calcula o operador BGK simétrico (modelo TRT)

void bgk_even ( double*, double*, double, double, double, double, double, double, double, double,
                double*, double*, double );

// Calcula o operador BGK antisimétrico (modelo TRT)

void bgk_odd ( double*, double*, double, double, double, double, double, double, double, double,
               double*, double*, double );

// Etapa de colisão usanto TRT

void trt_collision ( double*, double*, double, double, double, double, double, double*, double );

// Define os sítios para a etapa de propagação

void def_dir_prop ( double*, int*, int*, int, int, int );

// Etapa de propagação para a distribuição em um único sítio

void propag_part_site ( int*, double*, double*, int );

// Grava o campo de densidade

void rec_density ( int*, double*, double*, unsigned int, int, int, int, string );

// Grava o campo de velocidades (uma função distribuição)

void rec_velocity ( int*, double*, double*, unsigned int, int, int, int );

// Grava arquivo de recuperação

void rec_recovery ( double*, double*, int, int, int );

// Condição de contorno de derivada nula (Neumann)  direção x

void bound_dvnull_rho_x ( int*, double*, double*, int, double, int, int, int, double*, double );

//	Calcula a permeabilidade intrínseca

double conductivity ( double, double, double, double, double, int, int, int, int  );

// Retorna o momento na direção x

double quant_mov_x ( double*, double* );

// Calcula a densidade no sítio

double mass ( double* );

// Calcula um coeficiente angular (melhor reta)

double coefi ( double*, int, int, int, double* );

//===================== Inicio do programa principal =============================================//

int main ()
{

    //================= Inicialização de variáveis ===============================================//
    
    time_t tempo;
    
    int t_0 = ( int ) time( &tempo ); // tempo inicial em segundos

    string nome_geo;  // Nome do arquivo de geometria

    double ftesc;       // Dimensão do pixel

    int npassos;        // Número de passos de simulação

    int numarq;         // Número de arquivos a serem gravados

    double tau;         // Tempo de relaxação

    double visc;        // Viscosidade

    double rho_ini;      // Densidade (partículas/sítio)

	int max_threads = omp_get_max_threads();
	
	cout << "\nmax_threads = " << max_threads << endl;
	
	int set_threads = 0;
	
	cout << "\nNumber of threads: ";
	
	cin >> set_threads;
	
	omp_set_num_threads( set_threads );
        
    //----------------- Le o arquivo de inicialização --------------------------------------------//

    read_data ( nome_geo, ftesc, npassos, numarq, tau, visc, rho_ini );

    //--------------------------------------------------------------------------------------------//

    int intervalo;
    
    if ( numarq ) intervalo = npassos / numarq;
    
    else intervalo = npassos;

    int mais_passos = 0;

    unsigned int inicio = 0;

    unsigned int fim = npassos;

    //----------------- Específico para o modelo TRT ---------------------------------------------//

    double tau_sim = tau;

    double inv_tau = 1.0 / tau;

    double tau_ant = ( ( 8.0 - inv_tau ) / ( 8.0 * (  2.0 - inv_tau ) ) );

    //================= Le o arquivo de geometria do meio ========================================//

    int nx, ny, nz;

    ifstream fmatriz( nome_geo );

    string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    fmatriz.close();

    int *ini_meio = new int[ nx * ny * nz ];
    
    int ptos_meio = read_geo ( nome_geo, ini_meio, 0, 0 );

    double phi = ( double ) ptos_meio / ( double )( nx * ny * nz );

    cout << "\nPorosidade = " << phi << "    " << endl;

    //============ Aloca a memoria usada pelas outras matrizes ===================================//

    int tam_alloc = ptos_meio * nvel;

    double *ini_N = new double[tam_alloc];

    double *ini_N_novo = new double[tam_alloc];

    int *ini_dir = new int[tam_alloc];

    //============ Define Os Vetores de rede c_i =================================================//

    double *ini_c = new double[nvel*dim];

    def_lattice_d3q15 ( ini_c, W );

    //============ Define as direções de propagação ==============================================//

    def_dir_prop ( ini_c, ini_meio, ini_dir, nx, ny, nz );

    //================= Inicialização da rede ====================================================//

    double vx_ini = 0.0;

    double vy_ini = 0.0;

    double vz_ini = 0.0;

    char c;

    int c_int;

    do
    {
        cout << "\nRecuperar dados? (s/n): ";
        c = getchar();
        cout << "\r                        ";
        c_int = ( int ) c;
    }
    while ( c_int != 115 && c != 110 );

    //----------------- Inicialização padrão -----------------------------------------------------//

    if ( c_int == 110 )
    {
        for ( int pto = 0; pto < ptos_meio; pto++ )
        {

            double *N = ini_N + ( pto ) * nvel;

            dist_eq ( vx_ini, vy_ini, vz_ini, rho_ini, ini_c, N, W, one_over_c_s2 );
        }

        cout << endl;
    }

    //----------------- Inicialização com dados recuperados --------------------------------------//

    if ( c_int == 115 )
    {
        ifstream f_read ( "Arq_rec.dat" );

        f_read >> inicio;

        for ( int pto = 0; pto < ptos_meio; pto++ )
        {
            double vx, vy, vz, rho;

            f_read >> vx;

            f_read >> vy;

            f_read >> vz;

            f_read >> rho;

            double *N = ini_N + ( pto ) * nvel;

            dist_eq ( vx, vy, vz, rho, ini_c, N, W, one_over_c_s2 );
        }

        f_read.close();

        cout << endl;
    }

    //=================== Variáveis para o cálculo da permeabilidade =============================//
    
    double acc_x = 0.0;
    
    double acc_y = 0.0;

    double acc_z = 0.0;

    double soma_mx = 0.0;
    
    double delta_rho = 0.00001 * nx * rho_ini;
    
    double rho_in = rho_ini + delta_rho;
        
    double rho_out = rho_ini - delta_rho;
    
    double cs_2 = 1.0 / one_over_c_s2;
    
    double gradpress = - cs_2 * ( rho_out - rho_in ) / ( double ) nx;
		
    const double D_caract = ( double ) nx / 10.;
    
    double k_darcy = 0.0;
    
    double k_old = 0.0;
    
    double k_new = 0.0;
    
    unsigned int ptos_avg = 10;
    
    double* arq_k = new double[ ptos_avg ];
            
    double dk_dx = 1.0;
    
    double k_standart = 1.0;
    
    int contador = 0;
    
    //--------------------------------------------------------------------------------------------//

    unsigned int n_maior = ( int ) sqrt ( nx * nx + ny * ny + nz * nz );

    npassos = 10000 * n_maior;

    fim = npassos;

    //==================== Looping principal =====================================================//

loop:

    for (unsigned int passo = inicio; passo < fim; passo++ )
    {
        if ( passo % 2 == 0 ) cout << "\rStep : " << passo;

        //--------------- Imposição de condições de contorno -------------------------------------//

		#pragma omp parallel

		#pragma omp sections
        {
			#pragma omp section
			{
				int pos_x = 0;
									  
				bound_dvnull_rho_x ( ini_meio , ini_N, ini_c, pos_x, rho_in, nx, ny, nz, W,
									one_over_c_s2 );
			}

			#pragma omp section
			{
				int pos_x = nx - 1;

				bound_dvnull_rho_x ( ini_meio, ini_N, ini_c, pos_x, rho_out, nx, ny, nz, W,
									  one_over_c_s2 );
			}
		}

        //----------------------------------------------------------------------------------------//

		#pragma omp parallel for reduction(+:soma_mx)

        for ( int pto = 0; pto < ptos_meio; pto++ )
        {		
            //--------------- Aponta os ponteiros ------------------------------------------------//

            double *N = ini_N + ( pto ) * nvel;
            
            //--------------- Soma mx p/ calculo de k --------------------------------------------//
            
			//soma_mx = soma_mx + 0.5 * quant_mov_x ( N, ini_c );
			
            //--------------- Etapa de colisão ---------------------------------------------------//

            trt_collision ( N, ini_c, tau_sim, tau_ant, acc_x, acc_y, acc_z, W, one_over_c_s2 );

			//--------------- Soma mx p/ calculo de k --------------------------------------------//
            
			//soma_mx = soma_mx + 0.5 * quant_mov_x ( N, ini_c );
			
			soma_mx = soma_mx + quant_mov_x ( N, ini_c );
			
            //-------------- Etapa de propagacao -------------------------------------------------//

            propag_part_site ( ini_dir, ini_N, ini_N_novo, pto );
			
            //------------------------------------------------------------------------------------//
        }

        //-------------- Calcula e grava a permeabilidade ----------------------------------------//

        if ( passo % ( n_maior ) == 0 && passo > 0 )
        {
			
			double sum_grad = gradpress * ptos_meio * n_maior;
			
			k_darcy = conductivity ( soma_mx, sum_grad, ftesc, phi, visc, ptos_meio, nx, 
												passo, D_caract );

            //----------- Reinicializa variaveis -------------------------------------------------//

            soma_mx = 0.0;
            
            //------------------------------------------------------------------------------------//
        }

        //-------------- Grava os campos de velocidade e pressao ---------------------------------//

        if ( contador > 0 && contador == intervalo )
        {
			#pragma omp parallel

			#pragma omp sections
            {

				#pragma omp section
                {
                    rec_density ( ini_meio, ini_N, ini_c, passo, nx, ny, nz, "rho_" );
                }

				#pragma omp section
                {
                    rec_velocity ( ini_meio, ini_N, ini_c, passo, nx, ny, nz );
                }
            }

            contador = 0;
        }

        contador++;

        //------------------ Atualiza ( fnovo => f ) ---------------------------------------------//

        double *temp = ini_N;

        ini_N = ini_N_novo;

        ini_N_novo = temp;
        
        //-------------------- Verifica se deve terminar o cálculo -------------------------------//
        
        if ( passo % ( n_maior ) == 0 && passo > 0 )
        {
			if ( passo == n_maior ) k_standart = k_darcy;
			
			for ( int i = ( ptos_avg - 2 ); i >= 0; i-- ) arq_k[ i + 1 ] = arq_k[ i ];
			
			arq_k[ 0 ] = k_darcy;
			
			dk_dx = coefi ( arq_k, 0, ptos_avg, k_standart, &k_new );
		
			dk_dx = fabs( dk_dx );
		}
				        
        if ( passo % ( n_maior ) == 0 && passo > ( ptos_avg * n_maior ) )
        {
			double erro = 1.0;
			
			if ( k_new > 0.0 )
			{
				erro = fabs ( k_new - k_old ) / k_new;
			}
			else erro = fabs ( k_new - k_old );
			        
			if ( erro < 0.001 && dk_dx < 0.001 )
			{
				break;
			} 	
			
			k_old = k_new;		
		}

        //------------------ Grava arquivo de recuperação ----------------------------------------//

		int t_sec = ( int ) time( &tempo ) - t_0; 
    
        if ( t_sec > 3600 )
        {
            rec_recovery ( ini_N, ini_c, ptos_meio, passo, 8 );
            
            t_0 = ( int ) time( &tempo );
        }

    } //===================== Fim do looping principal ===========================================//

    cout << "\n\nMais quantos passos?" << endl << endl;
    cin >> mais_passos;

    if ( mais_passos != 0 )
    {
        inicio = fim;
        fim = inicio + mais_passos;
        goto loop;
    }
}

//================================================================================================//





//===================== Lê o arquivo de inicialização ============================================//
//
//      Input:
//      Output: nome do arquivo de geometria, tam. do pixel, passos, files, tempo de relaxação,
//              densidade inicial
//
//================================================================================================//

void read_data ( string& nome_geo, double& ftesc, int& npassos, int& numarq, double& tau,
                 double& visc, double& rho_ini )

{
    //--------------------------------------------------------------------------------------------//

    string nome_in = "data_in.txt";

    ifstream f_in( nome_in );

    //--------------------------------------------------------------------------------------------//

    string nome_out = "dat_out.txt";

    cout << "\nNome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    ofstream fdat( nome_out );

    fdat << "Nome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> nome_geo;

    cout << "\nNome do arquivo de geometria: " << nome_geo << endl;

    fdat << "Nome do arquivo de geometria: " << nome_geo << endl;

    //--------------------------------------------------------------------------------------------//

    string st_ft;
    
    f_in >> st_ft;

    f_in >> ftesc;  // Le a dimensao do dimensão do pixel ( m )

    cout << "\nDimensao do pixel = " << ftesc << " m" << endl;
 
    fdat << "\nDimensao do pixel = " << ftesc << " m" << endl;

    //--------------------------------------------------------------------------------------------//

    string st_steps;
    
    f_in >> st_steps;

    f_in >> npassos;

    cout << "\nNumero de passos: " << npassos << endl;

    fdat << "Numero de passos: " << npassos << endl;

    //--------------------------------------------------------------------------------------------//

    string st_files;
    
    f_in >> st_files;

    f_in >> numarq;

    cout << "\nNumero de arquivos: " << numarq << endl;

    fdat << "Numero de arquivos: " << numarq << endl;

    //--------------------------------------------------------------------------------------------//

    string st_tau;
    
    f_in >> st_tau;

    f_in >> tau;

    cout << "\nTempo de relaxacao: " << tau << endl;

    fdat << "Tempo de relaxacao: " << tau << endl;

    visc = 1.0 / 3.0 * ( tau - 0.5 );

    fdat << "viscosidade = " << visc << endl;

    //--------------------------------------------------------------------------------------------//

    string st_ro;
    
    f_in >> st_ro;

    f_in >> rho_ini;

    cout << "\nDensidade: " << rho_ini << endl;

    fdat << "Densidade: " << rho_ini << endl;

    //--------------------------------------------------------------------------------------------//

    f_in.close();

    fdat.close();
}

//================================================================================================//




//===================== Lê o arquivo de geometria ================================================//
//
//      Input: nome do arquivo de geometria, ponteiro para o meio
//      Output: número de pontos lidos
//
//================================================================================================//

int read_geo ( string nome_geo, int *meio, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nTamanho:  x = " << nx << ";  y = "   << ny << ";  z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//================================================================================================//




//===================== Calcula o produto interno ================================================//
//
//      Input: two vetors
//      Output: internal product
//
//================================================================================================//

double  prod_int ( double  vet_1[dim], double  vet_2[dim] )
{

    double  result = 0.0;

    for ( int i = 0; i < dim; i++ )
    {
        result = result + vet_1[i] * vet_2[i];
    }

    return result;
}

//================================================================================================//




//===================== Arredonda um número  =====================================================//
//
//      Input: double
//      Retorna: int
//
//================================================================================================//

int  sum_round ( double  num )
{
    int num_int;

    if ( num < 0 ) num_int = ( int ) ( num - 0.5 );
    else num_int = ( int ) ( num + 0.5 );

    return num_int;
}

//================================================================================================//




//===================== Define os vetores da rede D3Q15 ==========================================//
//
//      Input: pointer to the vectors
//      Output:
//
//================================================================================================//

void def_lattice_d3q15 ( double  *ini_c, double W[nvel] )
{
    double  *c;

    //------------- |c| = 0 ----------------------//

    c = ini_c;
    c[0] = 0;
    c[1] =  0;
    c[2] = 0;

    //------------- |c| = 1 ----------------------//

    c = ini_c + 1 * dim;
    c[0] =  1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 2 * dim;
    c[0] = -1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 3 * dim;
    c[0] =  0;
    c[1] =  1;
    c[2] =  0;

    c = ini_c + 4 * dim;
    c[0] =  0;
    c[1] = -1;
    c[2] =  0;

    c = ini_c + 5 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] =  1;

    c = ini_c + 6 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] = -1;

    //------------ |c| = sqrt(3) -----------------//

    c = ini_c + 7 * dim;
    c[0] =  1;
    c[1] =  1;
    c[2] =  1;

    c = ini_c + 8 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 9 * dim;
    c[0] = -1;
    c[1] =  1;
    c[2] =  1;

    c = ini_c + 10 * dim;
    c[0] =  1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 11 * dim;
    c[0] =  1;
    c[1] = -1;
    c[2] =  1;

    c = ini_c + 12 * dim;
    c[0] = -1;
    c[1] =  1;
    c[2] = -1;

    c = ini_c + 13 * dim;
    c[0] =  1;
    c[1] =  1;
    c[2] = -1;

    c = ini_c + 14 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] =  1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    W[0] = 2. / 9.;

    for ( int i = 1; i < 7; i++ ) W[i] = 1. / 9.;

    for ( int i = 7; i < nvel; i++ ) W[i] = 1. / 72.;
}

//================================================================================================//




//===================== Calcula a distribuição de equilíbrio =====================================//
//
//      Input: velocities, density, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//================================================================================================//

void dist_eq ( double vx, double vy, double vz, double rho, double *ini_c, double feq[nvel],
               double W[nvel], double one_over_c_s2 )
{
    double *c;

    double v[3];

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;

    double vquad = ( vx * vx + vy * vy + vz * vz );

    double cv;

    feq[0] = W[0] * rho * ( 1.0 - 0.5 * vquad * one_over_c_s2 );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        cv = prod_int ( c, v );

        feq[i] = W[i] * rho * ( 1.0 + cv * one_over_c_s2
                                + 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
                                - 0.5 * vquad * one_over_c_s2 );
    }
}

//================================================================================================//




//===================== Calcula densidade e velocidades ==========================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//================================================================================================//

void calcula ( double *f, double *ini_c, double& vx, double& vy, double& vz, double& rho )
{
    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;

    double one_over_rho;

    rho = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        double *c = ini_c + i * dim;

        rho = rho + f[i];

        mx = mx + c[0] * f[i];
        my = my + c[1] * f[i];
        mz = mz + c[2] * f[i];
    }

    if ( rho )
    {
        one_over_rho = 1.0 / ( rho );

        vx = mx * one_over_rho;
        vy = my * one_over_rho;
        vz = mz * one_over_rho;
    }
    else
    {
        vx = 0;
        vy = 0;
        vz = 0;
    }
}

//================================================================================================//




//===================== Calcula o operador BGK simétrico (modelo TRT) ============================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_even ( double f[nvel], double *ini_c, double vx, double vy, double vz, double rho,
                double acc_x, double acc_y, double acc_z, double tau, double op_col[nvel],
                double W[nvel], double one_over_c_s2 )
{
    double f_eq[nvel];

    double one_over_tau = 1.0 / tau;

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;

    dist_eq ( vx_acc, vy_acc, vz_acc, rho, ini_c, f_eq, W, one_over_c_s2 );

    op_col[0] = ( f_eq[0] - f[0] ) * one_over_tau;

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i+1] ) * 0.5;

        double f_eq_mais = ( f_eq[i] + f_eq[i+1] ) * 0.5;

        op_col[i] = ( f_eq_mais - f_mais ) * one_over_tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i-1] ) * 0.5;

        double f_eq_mais = ( f_eq[i] + f_eq[i-1] ) * 0.5;

        op_col[i] = ( f_eq_mais - f_mais ) * one_over_tau;
    }

}

//================================================================================================//




//===================== Calcula o operador BGK antisimétrico (modelo TRT) ========================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_odd ( double f[nvel], double *ini_c, double vx, double vy,
               double vz, double rho, double acc_x, double acc_y, double acc_z,
               double tau, double op_col[nvel], double W[nvel], double one_over_c_s2 )
{
    double f_eq[nvel];

    double one_over_tau = 1.0 / tau;

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;

    dist_eq ( vx_acc, vy_acc, vz_acc, rho, ini_c, f_eq, W, one_over_c_s2 );

    op_col[0] =  0.0;

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] - f[i+1] ) * 0.5;

        double f_eq_mais = ( f_eq[i] - f_eq[i+1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) * one_over_tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais =  ( f[i] - f[i-1] ) * 0.5;

        double f_eq_mais = ( f_eq[i] - f_eq[i-1] ) * 0.5;

        op_col[i] = ( f_eq_mais - f_mais ) * one_over_tau;
    }

}

//================================================================================================//




//===================== Etapa de colisão usando TRT ==============================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//================================================================================================//

void trt_collision ( double f[nvel], double *ini_c, double tau_sim, double tau_ant, double acc_x,
                     double acc_y, double acc_z, double W[nvel], double one_over_c_s2 )
{
    double vx, vy, vz, rho;

    calcula ( f, ini_c, vx, vy, vz, rho );
    
    //----------- Insere perdas de energia por efeito joule --------------------------//
    
    double f_lost = 0.125001; // regula a condutividade do material (água do mar = 0.97061)
    
    vx = f_lost * vx;
    vy = f_lost * vy;
    vz = f_lost * vz;

	//---------------------------------------------------------------------------------//
	
    double op_col_sim[nvel];

    bgk_even ( f, ini_c, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau_sim, op_col_sim, W,
               one_over_c_s2 );

    double op_col_ant[nvel];

    bgk_odd ( f, ini_c, vx, vy, vz, rho, acc_x, acc_y, acc_z, tau_ant, op_col_ant, W,
              one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col_sim[i] + op_col_ant[i];

    }
}

//================================================================================================//




//===================== Define os sítios para a etapa de propagação ==============================//
//
//      Input: lattice vectors, geometry, dimensions
//      Output: addresses, *ini_dir, and momentum lost, *ini_qlost
//
//================================================================================================//

void def_dir_prop ( double *ini_c, int *ini_meio, int *ini_dir, int nx, int ny, int nz )
{
    int i, j;
    double  *c;
    int *meio_local;
    int *meio_prop;

    int *dir;

    //---------------- Define os passos para propagação ------------------------------------------//

    double  step_x[nvel];
    double  step_y[nvel];
    double  step_z[nvel];

    int steps[nvel];

    for ( i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double  cx = c[0];
        double  cy = c[1];
        double  cz = c[2];

        int cx_int =  sum_round ( cx );
        int cy_int =  sum_round ( cy );
        int cz_int =  sum_round ( cz );

        steps[i] = abs ( cx_int );

        if ( abs ( cy_int ) > steps[i] )
        {
            steps[i] = abs ( cy_int );
        }

        if ( abs ( cz_int ) > steps[i] )
        {
            steps[i] = abs ( cz_int );
        }

        step_x[i] = cx / steps[i];
        step_y[i] = cy / steps[i];
        step_z[i] = cz / steps[i];
    }

    //-------------- Encontra as direções contrárias ---------------------------------------------//

    int i_op[nvel];

    for ( i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int cx_i = sum_round ( c[0] );
        int cy_i = sum_round ( c[1] );
        int cz_i = sum_round ( c[2] );

        for ( j = 1; j < nvel; j++ )
        {
            c = ini_c + j * dim;

            int cx_j = sum_round ( c[0] );
            int cy_j = sum_round ( c[1] );
            int cz_j = sum_round ( c[2] );

            if ( cx_j == -cx_i  && cy_j == -cy_i && cz_j == -cz_i )
            {
                i_op[i] = j;
            }
        }
    }

    //--------------------------------------------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio_local = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio_local )
                {
                    dir = ini_dir + ( *meio_local - 1 ) * nvel;

                    meio_prop = meio_local;

                    dir[0] = ( *meio_prop - 1 ) * nvel;

                    //----------------------------------------------------------------------------//

                    for ( i = 1; i < nvel; i++ )
                    {
                        int inv = 0; // número de inversões

                        double  stpx = step_x[i];
                        double  stpy = step_y[i];
                        double  stpz = step_z[i];

                        int x_prop = x;
                        int y_prop = y;
                        int z_prop = z;

                        double  x_prop_f = ( double ) x;
                        double  y_prop_f = ( double ) y;
                        double  z_prop_f = ( double ) z;

                        for ( int stp = 1; stp < steps[i] + 1; stp++ )
                        {
                            x_prop_f = x_prop_f + stpx;
                            y_prop_f = y_prop_f + stpy;
                            z_prop_f = z_prop_f + stpz;

                            x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                            y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                            z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                            meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                            if ( *meio_prop == 0 )
                            {
                                stpx = -stpx;
                                stpy = -stpy;
                                stpz = -stpz;

                                x_prop_f = x_prop_f + stpx;
                                y_prop_f = y_prop_f + stpy;
                                z_prop_f = z_prop_f + stpz;

                                x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                                y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                                z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                                meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                                inv++;
                            }
                        }

                        if ( inv % 2 == 0 )
                        {
                            dir[i] = i + ( *meio_prop - 1 ) * nvel;
                        }
                        else
                        {
                            dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
                        }
                    }

                    //----------------------------------------------------------------------------//
                }
            }
        }
    }
}

//================================================================================================//




//===================== Etapa de propagação para a um único sítio ================================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), momentum lost (*ini_qlost), point
//      Output: updated distribution function (*ini_f_new), summation of lost momentum
//
//================================================================================================//

void propag_part_site ( int *ini_dir, double *ini_f, double *ini_f_new, int pto )
{

    //--------------- Aponta os ponteiros --------------------------------------------------------//

    double *f = ini_f + ( pto ) * nvel;

    int *dir = ini_dir + ( pto ) * nvel;

    //--------------------------------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )
    {
        double *f_new;

        f_new = ini_f_new + dir[i];

        *f_new = f[i];
    }
}

//================================================================================================//




//===================== Grava o campo de velocidade - monofásico =================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_velocity ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                    int ny, int nz )
{
    char nomevel[50];

    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *f = ini_f + ( *meio - 1 ) * nvel;

                    double vx, vy, vz, rho;

                    calcula ( f, ini_c, vx, vy, vz, rho );

                    fvel << vx  << " " << vy << " " << vz << " ";
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//================================================================================================//




//===================== Grava o campo de densidade - um fluido ===================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_density ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                   int ny, int nz, string nomerho )
{
    nomerho = nomerho + to_string( passo ) + ".vtk";

    ofstream frho ( nomerho );

    frho << "# vtk DataFile Version 2.0" << endl;
    frho << "Densidade" << endl;
    frho << "ASCII" << endl;
    frho << "DATASET STRUCTURED_POINTS" << endl;
    frho << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    frho << "ASPECT_RATIO 1 1 1" << endl;
    frho << "ORIGIN 0 0 0" << endl;
    frho << "POINT_DATA " << nx*ny*nz << endl;
    frho << "SCALARS densidade double" << endl;
    frho << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *f = ini_f + ( *meio - 1 ) * nvel;

                    double rho = mass( f );

                    frho << rho << " ";
                }
                else
                {
                    frho << 0.0 << " ";
                }
            }

            frho << endl;
        }

    }
    frho.close();
}
//================================================================================================//



//===================== Grava arquivo de recuperação de dados ====================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_recovery ( double *ini_f, double *ini_c, int ptos_meio, int passo, int precision )
{
    cout << "\n\nGravando arquivo de recuperação..." << endl;

    ofstream f_rec ( "Arq_rec.dat", ios::binary );

    f_rec << passo << endl;

    for ( int pto = 0; pto < ptos_meio; pto++ )
    {
        double *f = ini_f + ( pto ) * nvel;

        double vx, vy, vz, rho;

        calcula ( f, ini_c, vx, vy, vz, rho );

        f_rec << setprecision( precision ) << vx << " ";

        f_rec << setprecision( precision ) << vy << " ";

        f_rec << setprecision( precision ) << vz << " ";

        f_rec << setprecision( precision ) << rho << " ";

    }

    f_rec.close();

    cout << "... ... ... !" << endl << endl;
}

//================================================================================================//




//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//
                           
void bound_dvnull_rho_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int x = posx;

    int infx = -1;
    
    if ( posx < nx / 2 ) infx = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int* meio_pto = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = ini_f + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = ini_f + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = mass( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    dist_eq ( 0., 0., 0., rho, ini_c, f_pto, W, one_over_c_s2 );
                }
            }
        }
    }
}

//================================================================================================//




//===================== Calcula a permeabilidade intrínseca ======================================//
//
//      Input: Quantidade de movimento do fluido (soma_mx), força total sobre o fluido (Q_lost),
//             dimensão do pixel (ftesc), porosidade (phi), viscosidade, dimensão x
//      Output: Permeabilidade
//
//================================================================================================//

double conductivity (double soma_mx, double Q_lost, double ftesc, double phi, double visc,
                              int ptos_meio, int nx, int passo, int D_caract )
{

    double mx_med = soma_mx / (double)( nx );

    //------------------------------------------------------------------------------------//

    double Q_lost_med = Q_lost / (double)( nx );
    //double Q_lost_med = sum_force / ( double ) ( nx );

    //------------------------------------------------------------------------------------//

    //phi = 1.0; // Teste
    double sigma = (ftesc)*(ftesc) * phi * visc * mx_med / Q_lost_med;

    //------------------------------------------------------------------------------------//

    cout << "\rStep : " << passo << "    Sigma = " << sigma
         << "	(S/m)" << endl << endl;

    //------------------------------------------------------------------------------------//

    double Re = ( D_caract ) * ( mx_med /  ptos_meio ) / visc;

    cout << "Reynolds = " << Re << endl << endl;

    ofstream fkm ("Sigma.dat",ios::app);
    fkm << passo << " " << sigma << endl;
    fkm.close();
    
    return sigma;
}

//================================================================================================//




//===================== Retorna o momento na direção x ===========================================//
//
//      Input: distribution function f[nvel], lattice vectors
//      Output: momentum in the x direction
//
//================================================================================================//

double quant_mov_x ( double f[nvel], double *ini_c )
{
    double *c;

    double mx = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        mx = mx + c[0] * f[i];
    }

    return mx;

}

//================================================================================================//



//===================== Calcula a densidade no sítio =============================================//
//
//      Input: distribution function
//      Output: density
//
//================================================================================================//

double mass ( double f[nvel] )
{

    double rho = 0.0;

    for ( int i = 0 ; i < nvel; i++ ) rho = rho + f[i];

    return rho;
}

//================================================================================================//




//===================== Calcula um coeficiente angular (melhor reta) =============================//
//
//      Input: 
//      Output: 
//
//================================================================================================//

double coefi ( double* y, int number , int npts, int k_std, double* avg )
{

	double sum_y = 0.0;

	double sum_x = 0.0;

	for ( int x = number; x < number + npts; x++ )
	{	
		sum_y = sum_y + y[x] / k_std;
		
		sum_x = sum_x + ( double ) x;
	}
    
	double x_avg = sum_x / ( double ) npts;
    
    double y_avg = sum_y / ( double ) npts;
    
    double sum_xy_yavg = 0.0;
    
    double sum_xx_xavg = 0.0;
    
    for ( int x = number; x < number + npts; x++ )
	{
		sum_xy_yavg = sum_xy_yavg + x * ( y[x] / k_std - y_avg );
		
		sum_xx_xavg = sum_xx_xavg + x * ( x - x_avg );
    }
            
	double dy_dx = sum_xy_yavg / sum_xx_xavg;
	
	*avg = y_avg * k_std;

	return dy_dx;
        
}
//================================================================================================//
